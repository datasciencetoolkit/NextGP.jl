module prepMatVec

using StatsModels, MixedModels, CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables, LinearAlgebra

import StatsModels.terms
StatsModels.terms!(path::String) = path #path for data and map

include("misc.jl")

export prep


"""
	function prep(f::StatsModels.TermOrTerms, inputData::DataFrame;userHints::Dict,path2ped,priorVCV)

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
function prep(f::StatsModels.TermOrTerms, inputData::DataFrame;userHints::Dict,path2ped,priorVCV)
	
#	any(typeof.(terms(f)).==ConstantTerm{Int64}) == false ? throw(ErrorException("Models without constant term are not allowed")) : nothing 
	
	terms4StatsModels = getTerms(f)

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
	iGRel = OrderedDict{Any,Any}()

	ME = Dict{Any,Any}()
	#ME = OrderedDict{Any,Any}()
	#map = OrderedDict{Any,Any}() 

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
		userData4ran = deepcopy(userData)
		userData4ran[!,[:ID,:Sire,:Dam]] .= pedigree[[findall(pedigree.origID.==x)[] for x in userData4ran.origID],[:ID,:Sire,:Dam]]
		
	end	

	#original id within pedigree
	#seemed to be IDs for only phenotyped ones????? from the ranMat()
		
	idRE = OrderedDict{Any,Any}()

	idFE = OrderedDict{Any,Any}() #fixed effects and their levels	

	#summarize input
	summarize = DataFrame(Variable=Any[],Term=Any[],Type=Any[],Levels=Int32[])

		

        for i in 1:length(f.rhs)
		if (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "SNP")
			(arg1,arg2,arg3...) = f.rhs[i].args_parsed
			arg1 = Symbol(repr(arg1))
			thisM = CSV.read(arg2,CSV.Tables.matrix,header=false)
			
			#str field can only be in GBLUP for marker related analysis
			if haskey(priorVCV,arg1) && in(:str,fieldnames(typeof(priorVCV[arg1])))
				RE[arg1] = Matrix(1.0*I,size(thisM,1),size(thisM,1))
				iGRel[arg1] = inv(makeG(thisM;method=priorVCV[arg1].type))
				push!(summarize,[arg1,"GBLUP",typeof(RE[arg1]),size(RE[arg1],2)])
			else

				thisM .-= mean(thisM,dims=1) 
				ME[arg1] = thisM
                       		thisM = 0 #I can directly merge to dict above
	#			map[arg1] = arg3[1]
				push!(summarize,[arg1,"SNP",typeof(ME[arg1]),size(ME[arg1],2)])
				iGRel[arg1] = [] ###temp
			end

			ME[:arg1] = Dict(:data=>thisM,:map=>arg3[1],:method=>"BayesPR",:str=>iGRel,:dims=>size(thisM),:levels=[M$i for i in 1:size(thisM,2)]) 


                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "PED")
                        arg = Symbol(repr((f.rhs[i].args_parsed)[1]))
			IDs,thisZ = ranMat(arg, :ID, userData, pedigree)
			RE[arg] = thisZ
			thisZ = 0
			idRE[arg] = [pedigree[findall(i.==pedigree.ID),:origID][] for i in IDs]
			push!(summarize,[arg,"PED",typeof(RE[arg]),size(RE[arg],2)])
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "|")
                        my_sch = schema(userData, userHints) #work on userData and userHints
			
			f.rhs[i].args_parsed[1] == ConstantTerm{Int64}(1) ? my_ApplySch = apply_schema(f.rhs[i].args_parsed[2], my_sch, MixedModels.MixedModel) : my_ApplySch = apply_schema(f.rhs[i], my_sch, MixedModels.MixedModel) 	
			#####ID is from the pheno  file directly, order not  checked!#####################################################
			arg1 = Symbol(repr((f.rhs[i].args_parsed)[1]))
                        arg2 = Symbol(repr((f.rhs[i].args_parsed)[2]))
			arg = Meta.parse(join([arg1,arg2]," | "))
                       	thisZ = modelcols(my_ApplySch, userData)
			RE[arg] = thisZ
			thisZ = 0
			idRE[arg] = unique(userData[!,arg2])
			push!(summarize,[arg,"|",typeof(RE[arg]),size(RE[arg],2)])

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

	ME = NamedTuple(ME)
	
	println("\n ---------------- Summary of input ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,nosubheader=true,alignment=:l)

	idFR = OrderedDict(:levelsFE => idFE, :levelsRE => idRE)


        return idFR, Ainv, iGRel, vec(yVec), FE, RE, ME
end

end

