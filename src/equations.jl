module equations

using StatsModels, MixedModels, CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables

include("misc.jl")

export mme

"""
        make_ran_matrix(x1::AbstractVector,x2::AbstractVector)

* Generates random effects matrix
* Initially works with onnly categorical vectors, to allow users add random effects as defined in StatsModels.jl
* So, the same varible can be used in two different function.
* For example, ran("dam","dam") can be similarly defined as (1|dam) for StatsModels.jl to create design matrices.

"""
function make_ran_matrix(x1::AbstractVector,x2::AbstractVector)
#        isa(x1, CategoricalArray) ||
#                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))
#        isa(x2, CategoricalArray) ||
#                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))

        u = sort(unique(x2));
        filter!(x->x≠0,u)
        Z = Matrix{Bool}(undef, length(x1), length(u))
        for i in eachindex(u)
        	@. Z[:, i] = x1 .== u[i]
        end
           return u,Z
       end


ranMat(arg1,arg2,data1,data2) = make_ran_matrix(data1[!,Symbol(arg1)],data2[!,Symbol(arg2)])


"""
	function mme(f::StatsModels.TermOrTerms, inputData::DataFrame;userHints::Dict,path2ped,paths2geno)

* `NextGP` relies on `StatsModels.jl` package for model expression (`f`), and fixed effect design matrix generation.
* Details for the model expression (`f`), and fixed effects coding specifications (e.g., effect or dummy coding) can be found at [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/).
* Design matrices for random effects are generated either own internal functions or using `StatsModels.jl`s `modelcols`, depending on how user defined the model term in the model.
* Reads in marker data, and mean-centers the columns.
* Finally returns lhs vector and rhs matrices.
* By default:
    * all `Int` rhs variables are made `Categorical`,
    * all `String` rhs variables (also those made `Categorical`) are dummy coded, except those defined by the user in `userHints`, 
    * all `Float` rhs variables are centered.
"""
function mme(f::StatsModels.TermOrTerms, inputData::DataFrame;userHints::Dict,path2ped,paths2geno)
        terms4StatsModels = String.(split(repr(f.rhs), ('+')))
        terms4StatsModels = replace.(terms4StatsModels, ":" => "")
        terms4StatsModels = [filter(x -> !isspace(x), trm) for trm in terms4StatsModels]

	#otherwise it changes original input data globally?????
	userData = deepcopy(inputData)

	for n in Symbol.(names(userData))
                if typeof(userData[!,n]).==Array{Int, 1}
                	userData[!,n] = CategoricalArray(userData[!,n])
        	end
        end


	for n in Symbol.(names(userData))
		if typeof(userData[!,n]).==Array{String, 1}
    			if !haskey(userHints,n)
				userHints[n] = StatsModels.DummyCoding()
			end
		end
	end


	#center cont. covariates	
	for n in Symbol.(names(userData))
		if n !== Symbol(repr(f.lhs))
        		if typeof(userData[!,n]).==Array{Float64, 1} || typeof(userData[!,n]).==Array{Float32, 1}
                		userData[!,n] .-= mean(userData[!,n],dims=1)
               		 end
		end
        end



        yVec = StatsModels.modelmatrix(f.lhs, userData)
	
        FE = OrderedDict{Any,Any}() #any to block work

        RE = OrderedDict{Any,Any}()

	ME = OrderedDict{Any,Any}()
	regionSizes = OrderedDict{String,Int64}()


        #read pedigree
	if isempty(path2ped)
		Ainv = []
	else

		pedigree,Ainv = makePed(path2ped,userData.ID)
		
		#sort data by pedigree. Needs to be carefully checked
		userData.origID = userData.ID
		userData.origSire = userData.Sire
		userData.origDam = userData.Dam
		userData.order = [findfirst(userData.origID .== x) for x in intersect(pedigree.origID,userData.origID)]
		sort!(userData, :order)
		select!(userData, Not(:order))
		
		#picking up new IDs (row/column number) from pedigree, and put into sire and dam in the phenotypic data 
		userData[!,[:ID,:Sire,:Dam]] .= pedigree[[findall(pedigree.origID.==x)[] for x in userData.origID],[:ID,:Sire,:Dam]]
		
	end	

	#original id within pedigree
	#seemed to be IDs for only phenotyped ones????? from the ranMat()
		
	idRE = OrderedDict{Any,Any}()

	idFE = OrderedDict{Any,Any}() #fixed effects and their levels	

	#summarize input
	summarize = DataFrame(Variable=Any[],Term=Any[],Type=Any[],Levels=Int32[])

        for i in 1:length(f.rhs)
		if (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "PR")
			arg1 = repr((f.rhs[i].args_parsed)[1])
			arg2 = parse(Int64,repr((f.rhs[i].args_parsed)[2]))
			path = paths2geno[Symbol(arg1)]
			thisM = CSV.read(path,CSV.Tables.matrix,header=false)
			thisM .-= mean(thisM,dims=1) 
			ME[arg1] = thisM
                        thisM = 0 #I can directly merge to dict above
			regionSizes[arg1] = arg2
			push!(summarize,[arg1,"BayesPR",typeof(ME[arg1]),size(ME[arg1],2)])
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "ran")
                        sym1 = repr((f.rhs[i].args_parsed)[1]) #now it is Symbol
                        sym2 = repr((f.rhs[i].args_parsed)[2]) #now it is from string			
			IDs,thisZ = ranMat(sym1, sym2, userData, pedigree)
			RE[(sym1,sym2)] = thisZ
			thisZ = 0
			idRE[(sym1,sym2)] = [pedigree[findall(i.==pedigree.ID),:origID][] for i in IDs]
			push!(summarize,[(sym1,sym2),"ran",typeof(RE[(sym1,sym2)]),size(RE[(sym1,sym2)],2)])
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "|")
                        my_sch = schema(userData, userHints) #work on userData and userHints
                        my_ApplySch = apply_schema(terms(f.rhs[i]), my_sch, MixedModels.MixedModel)
			#####NO IDs for this effect!!! Will be added later!!!!#####################################################
                       	thisZ = modelcols(my_ApplySch, userData)
			RE[terms4StatsModels[i]] = thisZ
			thisZ = 0
			push!(summarize,[f.rhs[i],"|",typeof(RE[terms4StatsModels[i]]),size(RE[terms4StatsModels[i]],2)])

                else
		my_sch = schema(userData, userHints)
		my_ApplySch = apply_schema(f.rhs[i], my_sch, MixedModels.MixedModel)
		idFE[terms4StatsModels[i]] = coefnames(my_ApplySch) 	
		thisX = modelcols(my_ApplySch, userData)
		FE[terms4StatsModels[i]] = thisX
		thisX = 0
		push!(summarize,[f.rhs[i],typeof(f.rhs[i]),typeof(FE[terms4StatsModels[i]]),size(FE[terms4StatsModels[i]],2)])
                end
        end
	
	println("\n ---------------- Summary of input ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,nosubheader=true,alignment=:l)

	idFR = OrderedDict(:levelsFE => idFE, :levelsRE => idRE)

        return idFR, Ainv, vec(yVec), FE, RE, ME, regionSizes
end

end

