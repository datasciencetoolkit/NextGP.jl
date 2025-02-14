module prepMatVec

using CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables, LinearAlgebra

#include("types.jl") already included in model.jl
include("model.jl")
include("misc.jl")
include("designMat.jl")

export prep


"""
	function prep(f::StatsModels.TermOrTerms, inputData::DataFrame;path2ped=[],priorVCV=[])

* `NextGP` relies on `StatsModels.jl` package for model expression (`f`), and fixed effect design matrix generation.
* Details for the model expression (`f`), and fixed effects coding specifications (e.g., effect or dummy coding) can be found at [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/).
* Design matrices for random effects are generated either own internal functions or using `StatsModels.jl`s `modelcols`, depending on how user defined the model term in the model.
* Reads in marker data, and mean-centers the columns.
* Finally returns lhs vector and rhs matrices.
* By default:
    * all `Int` rhs variables are made `Categorical`,
    * all `String` rhs variables (also those made `Categorical`) are dummy coded,
    * all `Float` rhs variables are centered.
"""
function prep(f, inputData::DataFrame;path2ped=[],priorVCV=[])
	
#	any(typeof.(terms(f)).==ConstantTerm{Int64}) == false ? throw(ErrorException("Models without constant term are not allowed")) : nothing 
	
	modelTerms = getTerms(f)

	userData = deepcopy(inputData)

	for n in Symbol.(names(userData))
                if typeof(userData[!,n]).==Array{Int, 1}
                	userData[!,n] = CategoricalArray(userData[!,n])
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


	#only for single trait now
        yVec = makeX(userData,f.lhs)[:data]
	
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

		

        for (k,v) in modelTerms
		println("$k is a $(typeof(v))")
		if isa(v,GenomicTerm)			
			thisM = CSV.read(v.path,CSV.Tables.matrix,header=false,delim=' ') #now white single white space is used 
			#drops cols if any value is missing. Later should check map files etc..
			thisM = thisM[:,.!(any.(ismissing, eachcol(thisM)))]
			#
			thisM = Matrix{Float64}(thisM)
			isempty(v.map) ? nowMap=[] : nowMap=v.map
			
			#str field can only be in GBLUP for marker related analysis
			if haskey(priorVCV,k) && in(:str,fieldnames(typeof(priorVCV[k])))
				iGRel = Symmetric(inv(makeG(thisM;method=priorVCV[k].type)))
				push!(summarize,[k,"GBLUP",typeof(iGRel),size(iGRel,2)])
				Z[k] = Dict(:data=>Matrix(1.0*I,size(thisM,1),size(thisM,1)),:map=>nowMap,:method=>"GBLUP",:str=>"G",:iVarStr=>iGRel,:dims=>size(iGRel),:levels=>["Ind$i" for i in 1:size(thisM,2)]) 	
		
			else
				thisM .-= mean(thisM,dims=1)
				M[k] = Dict(:data=>thisM,:map=>nowMap,:method=>"SNP",:str=>"I",:iVarStr=>[],:dims=>size(thisM),:levels=>["M$i" for i in 1:size(thisM,2)]) 			
				push!(summarize,[k,"Marker Effect",typeof(thisM),size(thisM,2)])
			end
			thisM = 0
		elseif isa(v,PedigreeTerm)
			IDs,thisZ = ranMat(k, :ID, userData4ran, pedigree)
			ids = [pedigree[findall(i.==pedigree.ID),:origID][] for i in IDs]
			Z[k] = Dict(:data=>thisZ,:method=>"BLUP",:str=>"A",:iVarStr=>Ainv,:dims=>size(Ainv),:levels=>ids) 	
			push!(summarize,[k,"PED",typeof(thisZ),size(thisZ,2)])
			thisZ = 0                
                else    
			#if isa(v,ConstantTerm)
			#	X[k] = Dict(:data=>ones(size(userData,1)),:map=>[],:method=>"FixedEffects",:nCol=>1,:levels=>"Intercept")
			#elseif isa(v,DataTerm)
			#	X[k] = makeX(userData,k)
			#elseif isa(v,FunctionTerm)
			#	X[k] = makeX(userData,modelTerms[k].cols)
			#	X[k][:data] = map(getproperty(Main, modelTerms[k].fun),X[k][:data])
			#elseif isa(v,InteractionTerm)
			#	X[k] = makeX(userData,modelTerms[k].cols)
			#else nothing
			#end
			
			X[k] = designMat(k,v,userData)
			push!(summarize,[k,typeof(k),typeof(X[k][:data]),X[k][:nCol]])
                end
        end

	
	println("\n ---------------- Summary of input ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,alignment=:l)

        return yVec, X, Z, M
end

end

