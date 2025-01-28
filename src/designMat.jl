
using DataFrames


#should work for only categorical variables
function makeXCat(tempData::Vector)
	#println(tempData)
	colLevels = unique(tempData)
	#println(colLevels)
  	dictCol = Dict()
	for (i,c) in enumerate(colLevels)
  		dictCol[c] = i
	end
	ii = 1:size(tempData,1)            # get row numbers
	jj = [dictCol[i] for i in tempData]  # get column numbers using list comprehension
	codedCol = Matrix(sparse(ii,jj,1.0))	
	return dictCol, codedCol
end


#Int variable is considered as categorical
#no Int String interaction yet
function makeX(df::DataFrame,col::Symbol)
	#println("processing $col")
	tempData = df[!,col]
	if isa(tempData,Vector{String})
		#println("tempData is a String")
		Xnow = makeXCat(tempData)
	elseif isa(tempData,Vector{Float64})
		#println("tempData is a Float")
		Xnow = tempData
	elseif isa(tempData,Vector{Int64})
		#println("tempData is a Int")
		Xnow = makeXCat(tempData)
	end
	return Xnow
end


#should work for only categorical variables
function makeXInt(tempData::Matrix)
	#println(tempData)
	colLevels = unique(tempData,dims=1)
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
	return dictCol, codedCol
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
		Xnow = prod.(eachrow(tempData))
	elseif isa(tempData,Matrix{Int64})
		#println("tempData is a MAT Int")
		Xnow = makeXInt(tempData)
	end
	return Xnow
end
