#module prepMatVec

using CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables, LinearAlgebra

include("misc.jl")
include("designMat.jl")

#all int are made categorical
function prepData!(inputData,formula)
	#make in categorical
	#for n in Symbol.(names(inputData))
	#	if isa(inputData[!,n],Vector{Int})
        #        	inputData[!,n] = CategoricalArray(inputData[!,n])
        #	end
        #end

	
	#center cont. covariates CURRENTLY ONLY APPLIES TO SINGLE-TRAIT ANALYSIS
	for n in Symbol.(names(inputData))
		if n !== formula.lhs
        		if typeof(inputData[!,n]).==Array{Float64, 1} || typeof(inputData[!,n]).==Array{Float32, 1}
                		inputData[!,n] .-= mean(inputData[!,n],dims=1)
               		 end
		end
        end
	return inputData
end

#can modify inputData #Need to avoid modifications as I read the ped multiple times given multiple random effects read ot multiple times
#no ID,Sire,Dam in the pedigree File
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
	return inputData,pedigree,Ainv
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
	#Do I need OrderedDict anywhere for Tuple() related matrix for E
	modelInformation = Dict{Any,Any}()
	X = OrderedDict{Any,Any}()
	Z = OrderedDict{Any,Any}()
	M = OrderedDict{Any,Any}()
	E = OrderedDict{Any,Any}()
	
	if length(f) == 1 #both traits have the same model terms
		modelRHSTerms = getRHSTerms(f[1])
		modelLHSTerms = getLHSTerms(f[1])
		#yVec is a vector if one response variable, matrix otherwise. functions.jl may need to be changed to work with matrix yCorr also.
		if length(modelLHSTerms) == 1
			inputData = CSV.read(f[1].data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
			inputData = prepData!(inputData,f[1])
			Y = makeX(inputData,f[1].lhs)[:data]
			E[f[1].lhs] = Dict{Any,Any}()
			modelInformation[collect(keys(modelLHSTerms))[]] = modelRHSTerms#keys(modelRHSTerms)
		elseif length(modelLHSTerms) > 1
			inputData = CSV.read(f[1].data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
			Y = hcat([makeX(inputData,k)[:data] for (k,v) in modelLHSTerms]...)
			E[Tuple(collect(keys(modelLHSTerms)))] = Dict{Any,Any}()
			modelInformation[Tuple(collect(keys(modelLHSTerms)))] = modelRHSTerms #keys(modelRHSTerms)
		end
	elseif length(f) > 1
		Y = [] #Matrix(undef,0,length(f))
		modelLHSTerms = Dict()
		modelRHSTerms = Dict()
		for (i,fi) in enumerate(f)
			println("reading $i $fi")
			inputData = CSV.read(fi.data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
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
	#println("inputData $inputData")
	#println(E)

	
	#summarize input
	summarize = DataFrame(Variable=Any[],Term=Any[],Type=Any[],Levels=Int32[])

	pedigree = []
	Ainv = []

        for (k,v) in modelRHSTerms
		if isa(v,GenomicTerm)			
			isempty(v.map) ? nowMap=[] : nowMap=v.map
			if haskey(priorVCV,k) && isa(priorVCV[k],GBLUPType) && isempty(priorVCV[k].G)
				println("GBLUP WITH COMPUTED G")
				thisM = CSV.read(String(v.path),CSV.Tables.matrix,header=false,delim=' ') #now white single white space is used 
				#drops cols if any value is missing. Later should check map files etc..
				thisM = thisM[:,.!(any.(ismissing, eachcol(thisM)))]
				#
				thisM = Union{Int64,Matrix{Float64}}(thisM)
				thisM = map(Int, thisM[:,1]) #I convert first column to Int to make sure it is correctly matches with IDs in the input data
				
				ids,thisZ = make_ran_matrix(inputData[!,:ID],thisM[:,1])
				G = makeG(thisM[:,2:end];method=priorVCV[k].methodG) #ignore my fake numerical IDs
				Z[k] = Dict(:data=>thisZ,:map=>nowMap,:method=>"GBLUP",:str=>"G",:iVarStr=>inv(G),:dims=>size(thisM,1),:levels=>ids) 	
				push!(summarize,[k,"GBLUP",typeof(Z[k][:iVarStr]),size(Z[k][:iVarStr],1)])
			elseif haskey(priorVCV,k) && isa(priorVCV[k],GBLUPType) && !isempty(priorVCV[k].G)
				println("GBLUP WITH GIVEN G")
				if isa(priorVCV[k].G,Matrix{Float64})
					ids,thisZ = make_ran_matrix(inputData[!,:ID],priorVCV[k].G[:,1]) #first column must be integer ID!! Currently not very flexible
					Z[k] = Dict(:data=>thisZ,:map=>nowMap,:method=>"GBLUP",:str=>"G",:iVarStr=>inv(priorVCV[k].G),:dims=>size(priorVCV[k].G),:levels=>ids)
				elseif isa(priorVCV[k].G,String) #not working because of type.jl part. GBLUP with string not exported for some reason. Fix later
					G = CSV.read(priorVCV[k].G,CSV.Tables.matrix,header=false,delim=' ') #now white single white space is used
					ids,thisZ = make_ran_matrix(inputData[!,:ID],G[:,1]) #first column must be integer ID!! Currently not very flexible
					Z[k] = Dict(:data=>thisZ,:map=>nowMap,:method=>"GBLUP",:str=>"G",:iVarStr=>inv(G),:dims=>size(G),:levels=>ids)
					G = 0
				end
				push!(summarize,[k,"GBLUP",typeof(Z[k][:iVarStr]),size(Z[k][:iVarStr],1)])
			elseif haskey(priorVCV,k) && isa(priorVCV[k],RandomMarkerEffect)
				println("MARKER EFFECT TYPE")
				thisM = CSV.read(String(v.path),CSV.Tables.matrix,header=false,delim=' ') #now white single white space is used 
				#drops cols if any value is missing. Later should check map files etc..
				thisM = thisM[:,.!(any.(ismissing, eachcol(thisM)))]
				#
				thisM = Matrix{Float64}(thisM)
				
				thisM .-= mean(thisM,dims=1)
				M[k] = Dict(:data=>thisM,:map=>nowMap,:method=>"SNP",:str=>"I",:iVarStr=>[],:dims=>size(thisM),:levels=>["M$i" for i in 1:size(thisM,2)]) 			
				push!(summarize,[k,"Marker Effect",typeof(thisM),size(thisM,2)])
			else throw(ArgumentError("Could not understand the analysis method for $k"))
			end
			thisM = 0
			thisZ = 0
		elseif isa(v,PedigreeTerm)
			if isempty(pedigree)
				println("pedigree is being computed for $k")
				inputData,pedigree,Ainv = usePedigree!(v.path,inputData)
				println("Ainv: $Ainv")
			else 
				println("pedigree is defined already for k")#nothing
			end
			IDs,thisZ = ranMat(k, :ID, inputData, pedigree)
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

#end

