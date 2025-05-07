module prepMatVec

using CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables, LinearAlgebra

include("model.jl")
include("misc.jl")
include("designMat.jl")
include("types.jl")

export @model,prep

function prepData!(inputData,f)
	#make in categorical
	for n in Symbol.(names(inputData))
		if isa(inputData[!,n],Vector{Int})
                	inputData[!,n] = CategoricalArray(inputData[!,n])
        	end
        end

	#center cont. covariates	
	for n in Symbol.(names(inputData))
		if n !== Symbol(repr(f.lhs))
        		if typeof(inputData[!,n]).==Array{Float64, 1} || typeof(inputData[!,n]).==Array{Float32, 1}
                		inputData[!,n] .-= mean(inputData[!,n],dims=1)
               		 end
		end
        end
	return inputData
end

#can modify inputData
function usePedigree!(path2ped,inputData)
	#read pedigree
	if isempty(path2ped)
		Ainv = []
	else

		pedigree,Ainv = makePed(path2ped,inputData.ID)
		Ainv = Symmetric(Ainv)
		
		#sort data by pedigree. Needs to be carefully checked
		inputData.origID = inputData.ID
		inputData.origSire = inputData.Sire
		inputData.origDam = inputData.Dam
		inputData.order = [findfirst(inputData.origID .== x) for x in intersect(pedigree.origID,inputData.origID)]
		sort!(inputData, :order)
		select!(inputData, Not(:order))
		
		#picking up new IDs (row/column number) from pedigree, and put into sire and dam in the phenotypic data
		userData4ran = deepcopy(inputData)
		inputData[!,[:ID,:Sire,:Dam]] .= pedigree[[findall(pedigree.origID.==x)[] for x in userData4ran.origID],[:ID,:Sire,:Dam]]
		
	end	

	#original id within pedigree
	#seemed to be IDs for only phenotyped ones????? from the ranMat()	
	#idRE = OrderedDict{Any,Any}()
	return inputData,Ainv
end


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
function prep(f;path2ped=[],priorVCV=[]) ### THE REST OF THE CODE FOR XZM SHOUld also come here, otherwise input data is only the last one in the memory!

	#I am assuming different variables names for each model. y1 = a1 + b1 +... and y2 = a2 + b2 +....
	modelInformation = Dict{Any,Any}()
	X = Dict{Any,Any}()
	Z = Dict{Any,Any}()
	M = Dict{Any,Any}()
	E = Dict{Any,Any}()
	
	if length(f) == 1 #both traits have the same model terms
		modelRHSTerms = getRHSTerms(f[1])
		modelLHSTerms = getLHSTerms(f[1])
		#yVec is a vector if one response variable, matrix otherwise. functions.jl may need to be changed to work with matrix yCorr also.
		if length(modelLHSTerms) == 1
			inputData = CSV.read(f[1].data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
			inputData = prepData!(inputData,f[1])
			inputData,Ainv = usePedigree!(path2ped,inputData)
			Y = makeX(inputData,f[1].lhs)[:data]
			E[f[1].lhs] = Dict{Any,Any}()
			modelInformation[collect(keys(modelLHSTerms))[]] = modelRHSTerms#keys(modelRHSTerms)
		elseif length(modelLHSTerms) > 1
			inputData = CSV.read(f[1].data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
			Y = hcat([makeX(inputData,k)[:data] for (k,v) in modelLHSTerms]...)
			#[E[k] = Dict{Any,Any}() for (k,v) in modelLHSTerms]
			E[Tuple(collect(keys(modelLHSTerms)))] = Dict{Any,Any}()
			modelInformation[collect(keys(modelLHSTerms))] = keys(modelRHSTerms)
		end
	elseif length(f) > 1
		Y = [] #Matrix(undef,0,length(f))
		modelLHSTerms = Dict()
		modelRHSTerms = Dict()
		for (i,fi) in enumerate(f)
			println("reading $i $fi")
			inputData = CSV.read(fi.data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
			inputData,Ainv = usePedigree!(path2ped,inputData)
			LHSfi = getLHSTerms(fi)
			println("LHSfi: $LHSfi")
			RHSfi = getRHSTerms(fi)
			modelLHSTerms = merge!(modelLHSTerms,LHSfi)
			modelRHSTerms = merge!(modelRHSTerms,RHSfi)
			yi = makeX(inputData,collect(keys(LHSfi))[])[:data]
			push!(Y,yi)
			E[collect(keys(LHSfi))[]] = Dict{Any,Any}()
			modelInformation[collect(keys(LHSfi))[]] = collect(keys(RHSfi))
		end
		Y = hcat(Y...)
	else throw(ArgumentError("model expression is not valid"))
	end

	println("modelInformation: $(modelInformation)")

	#println(modelLHSTerms)
	#println(modelRHSTerms)
	#println(inputData)
	println(E)

	
	#summarize input
	summarize = DataFrame(Variable=Any[],Term=Any[],Type=Any[],Levels=Int32[])

        for (k,v) in modelRHSTerms
		if isa(v,GenomicTerm)			
			thisM = CSV.read(String(v.path),CSV.Tables.matrix,header=false,delim=' ') #now white single white space is used 
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
			X[k] = designMat(k,v,inputData)
			push!(summarize,[k,typeof(k),typeof(X[k][:data]),X[k][:nCol]])
                end
        end

	
	println("\n ---------------- Summary of input ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,alignment=:l)

        return Y, X, Z, M, E, modelInformation
end

end

