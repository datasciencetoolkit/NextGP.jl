
using DataFrames, SparseArrays

#remove any column from data as reference column in dummy coding
dropCol(matrix::AbstractMatrix,j) = matrix[:, deleteat!(collect(axes(matrix, 2)), j)]
#ref level is always first now! Later reflevel shouild be changed to k!=reflevel with k being real level not coded value. Also drop col should be adapted sa such

function designMat(k,v,userData)
	if isa(v,ConstantTerm)
		println("$v is a ConstantTerm")
		X0 = Dict(:data=>ones(size(userData,1)),:map=>[],:method=>"FixedEffects",:nCol=>1,:levels=>"Intercept")
	elseif isa(v,DataTerm)
		println("$v is a DataTerm")
		X0 = makeX(userData,k)
	elseif isa(v,FunctionTerm)
		X0 = makeX(userData,modelTerms[k].cols)
		X0[:data] = map(getproperty(Main, modelTerms[k].fun),X[k][:data])
	elseif isa(v,InteractionTerm)
		X0 = makeX(userData,modelTerms[k].cols)
	else nothing
	end
	return X0
end

function dropLevel(levelCodesDict,matrix;refKey=[])
	refKey = isempty(refKey) ? first(keys(levelCodesDict)) : refKey
	refValue = levelCodesDict[refKey]
	levelCodesDict = sort(Dict(k => v for (k, v) in levelCodesDict if k != refKey))
	matrix = dropCol(matrix,refValue)
	return levelCodesDict,matrix 
end

#maybe not needed Union thing below?

#should work for only categorical variables
function makeXCat(tempData::Vector,col::Symbol)
	#println(tempData)
	#makes values and keys ordered as wanted
	colLevels = sort(unique(tempData))
  	dictCol = Dict()
	for (i,c) in enumerate(colLevels)
  		dictCol[c] = i
	end
	ii = 1:size(tempData,1)            # get row numbers
	jj = [dictCol[i] for i in tempData]  # get column numbers using list comprehension
	codedCol = Matrix(sparse(ii,jj,1.0))
	dictCol = sort(dictCol,byvalue=true)
	###DUMMY CODING
	dictCol,codedCol = dropLevel(dictCol,codedCol;refKey=[])
	return Dict(:data=>codedCol,:map=>[],:method=>"FixedEffects",:nCol=>size(codedCol,2),:levels=>"$col"*":".*keys(dictCol))
end


#Int variable is considered as categorical
#no Int String interaction yet
function makeX(df::DataFrame,col::Symbol)
	#println("processing $col")
	tempData = df[!,col]
	if isa(tempData,Vector{String})
		#println("tempData is a String")
		Xnow = makeXCat(tempData,col)
	elseif isa(tempData,Vector{Float64})
		#println("tempData is a Float")
		Xnow = Dict(:data=>tempData,:map=>[],:method=>"FixedEffects",:nCol=>1,:levels=>String(col))
	elseif isa(tempData,Vector{Int64})
		#println("tempData is a Int")
		Xnow = makeXCat(tempData)
	end
	return Xnow
end


#should work for only categorical variables
function makeXInt(tempData::Matrix)
	#println(tempData)
	colLevels = sort(unique(tempData,dims=1))
	colLevels = Array.(eachrow(colLevels))
	#println(colLevels)
  	dictCol = Dict()
	for (i,c) in enumerate(colLevels)
  		dictCol[c] = i
	end
	#println(dictCol)
	ii = 1:size(tempData,1)            # get row numbers
	jj = [dictCol[i] for i in eachrow(tempData)]  # get column numbers using list comprehension
	codedCol = Matrix(sparse(ii,jj,1.0))
	dictCol = sort(dictCol,byvalue=true)
	dictCol,codedCol = dropLevel(dictCol,codedCol;refKey=[])
	return Dict(:data=>codedCol,:map=>[],:method=>"FixedEffects",:nCol=>length(colLevels),:levels=>keys(dictCol))
end


#Int variable is considered as categorical
#no Int String interaction yet
function makeX(df::DataFrame,col::Vector{Symbol})
	#println("processing $col")
	tempData = Matrix(df[!,col])	
	if isa(tempData,Matrix{String})
		#println("tempData is a MAT String")
		Xnow = makeXInt(tempData)
	elseif isa(tempData,Matrix{Float64})
		#println("tempData is a MAT Float")
		Xnow = Dict(:data=>prod.(eachrow(tempData)),:map=>[],:method=>"FixedEffects",:nCol=>1,:levels=>String.(col))
	elseif isa(tempData,Matrix{Int64})
		#println("tempData is a MAT Int")
		Xnow = makeXInt(tempData)
	end
	return Xnow
end
