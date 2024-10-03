__precompile__(false) 

module prepMatVec

using StatsModels, MixedModels, CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables, LinearAlgebra

#import StatsModels.terms
#StatsModels.terms!(path::String) = path #path for data and map
import StatsModels.parse!
parse!(path::String, protected) = path
StatsModels.termvars(path::String) = path #path for data and map


include("misc.jl")

export prep


"""
	function prep(f::StatsModels.TermOrTerms, inputData::DataFrame;userHints=Dict{Symbol,Any}(),path2ped=[],priorVCV=[])

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
function prep(f::StatsModels.TermOrTerms, inputData::DataFrame;userHints=Dict{Symbol,Any}(),path2ped=[],priorVCV=[])
	
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
	
	X = Dict{Any,Any}()
	Z = Dict{Any,Any}()
	M = Dict{Any,Any}()

        #read pedigree
	if isempty(path2ped)
		Ainv = []
	else

		pedigree,Ainv = makePed(path2ped,userData.ID)
		Ainv = Symmetric(Ainv)
		
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

	#summarize input
	summarize = DataFrame(Variable=Any[],Term=Any[],Type=Any[],Levels=Int32[])

		

        for i in 1:length(f.rhs)
		if (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "SNP")
			(arg1,arg2,arg3...) = f.rhs[i].args_parsed
			arg1 = Symbol(repr(arg1))
			thisM = CSV.read(arg2,CSV.Tables.matrix,header=false,delim=' ') #now white single white space is used 
			#drops cols if any value is missing. Later should check map files etc..
			thisM = thisM[:,.!(any.(ismissing, eachcol(thisM)))]
			#
			thisM = Matrix{Float64}(thisM)
			
			#str field can only be in GBLUP for marker related analysis
			if haskey(priorVCV,arg1) && in(:str,fieldnames(typeof(priorVCV[arg1])))
				iGRel = Symmetric(inv(makeG(thisM;method=priorVCV[arg1].type)))
				push!(summarize,[arg1,"GBLUP",typeof(iGRel),size(iGRel,2)])
				Z[arg1] = Dict(:data=>Matrix(1.0*I,size(thisM,1),size(thisM,1)),:map=>arg3[1],:method=>"GBLUP",:str=>"G",:iVarStr=>iGRel,:dims=>size(iGRel),:levels=>["Ind$i" for i in 1:size(thisM,2)]) 	
		
			else
				thisM .-= mean(thisM,dims=1)
				isempty(arg3) ? nowMap=[] : nowMap=arg3[1]
				M[arg1] = Dict(:data=>thisM,:map=>nowMap,:method=>"SNP",:str=>"I",:iVarStr=>[],:dims=>size(thisM),:levels=>["M$i" for i in 1:size(thisM,2)]) 			
				push!(summarize,[arg1,"Marker Effect",typeof(thisM),size(thisM,2)])
			end
			thisM = 0

                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "PED")
                        arg = Symbol(repr((f.rhs[i].args_parsed)[1]))
			IDs,thisZ = ranMat(arg, :ID, userData4ran, pedigree)
			ids = [pedigree[findall(i.==pedigree.ID),:origID][] for i in IDs]
			Z[arg] = Dict(:data=>thisZ,:method=>"BLUP",:str=>"A",:iVarStr=>Ainv,:dims=>size(Ainv),:levels=>ids) 	
			push!(summarize,[arg,"PED",typeof(thisZ),size(thisZ,2)])
			thisZ = 0
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "|")
                        my_sch = schema(userData, userHints) #work on userData and userHints
			
			f.rhs[i].args_parsed[1] == ConstantTerm{Int64}(1) ? my_ApplySch = apply_schema(f.rhs[i].args_parsed[2], my_sch, MixedModels.MixedModel) : my_ApplySch = apply_schema(f.rhs[i], my_sch, MixedModels.MixedModel) 	
			#####ID is from the pheno  file directly, order not  checked!#####################################################
			arg1 = Symbol(repr((f.rhs[i].args_parsed)[1]))
                        arg2 = Symbol(repr((f.rhs[i].args_parsed)[2]))
			arg = Meta.parse(join([arg1,arg2]," | "))
                       	thisZ = modelcols(my_ApplySch, userData)
			ids = unique(userData[!,arg2])
			strI = Matrix(1.0*I,size(thisZ,2),size(thisZ,2))
			Z[arg] = Dict(:data=>thisZ,:method=>"|",:str=>"I",:iVarStr=>strI,:dims=>size(strI),:levels=>ids) 	
			push!(summarize,[arg,"|",typeof(thisZ),size(thisZ,2)])
			thisZ = 0
                else
			my_sch = schema(userData[!,intersect(Symbol.(names(userData)),terms4StatsModels)],userHints) #can be done only once above
#			my_sch = schema(userData[!,[Symbol(f.rhs[i])]]) #will crash for general mean
#			my_sch = schema(userData, userHints)
			my_ApplySch = apply_schema(f.rhs[i], my_sch, MixedModels.MixedModel)
			levelX = coefnames(my_ApplySch)
			thisX = modelcols(my_ApplySch, userData)
			isa(thisX,Vector) ? nCol = 1 : nCol = size(thisX,2)
			X[terms4StatsModels[i]] = Dict(:data=>thisX,:map=>[],:method=>"FixedEffects",:nCol=>nCol,:levels=>levelX) 
			push!(summarize,[f.rhs[i],typeof(f.rhs[i]),typeof(thisX),nCol])
			thisX = 0
                end
        end

	
	println("\n ---------------- Summary of input ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,alignment=:l)

        return vec(yVec), X, Z, M
end

end

