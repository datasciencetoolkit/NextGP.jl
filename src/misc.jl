using Printf
using DataFrames
using CSV
using PedigreeBase
using StatsBase


function getTerms(f)
        terms4StatsModels = String.(split(repr(f.rhs), ('+')))
        terms4StatsModels = replace.(terms4StatsModels, ":" => "")
        terms4StatsModels = [filter(x -> !isspace(x), trm) for trm in terms4StatsModels]
        terms4StatsModels = Symbol.(terms4StatsModels)
        return(terms4StatsModels)
end


"""
        make_ran_matrix(x1::AbstractVector,x2::AbstractVector)

* Generates random effects matrix
* Initially works with onnly categorical vectors, to allow users add random effects as defined in StatsModels.jl

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


ranMat(arg1,arg2,data1,data2) = make_ran_matrix(data1[!,arg1],data2[!,arg2])


"""
	 function(dataLHS::DataFrame)
Makes LHS of formula expression, when the LHS is a Matrix
It should be run by the user before @formula as `getLhs(myData)`
Then `f = @eval @formula(\$Y ~ 1 + x)`
"""
function getLhs(dataLHS::DataFrame)
	varNames = Symbol.(names(select(dataLHS, Not(:x))))
	expLHS = Expr(:call, :+, (varNames...))
	return expLHS
end

"""
Makes LHS of formula, when the LHS is a Matrix
It should be run by the user before @formula as `getLhs2(myData)`
Then `f = Y ~ @formula(0 ~ 1 + x).rhs`
"""
function getLhs2(dataLHS::DataFrame)
        dataLHS = select(dataLHS, Not(:x))
        Y = sum(term.(propertynames(dataLHS)))
        return Y
end


"""
	makeA(s::Any, d::Any)
Makes pedigree-based relationship matrix.
adapted from http://morotalab.org/Mrode2005/relmat/createA.txt

"""
function makeA(s::Any, d::Any)
    s = convert(Vector{Int64},s)
    d = convert(Vector{Int64},d)
    n = length(s)
    N = n + 1
    A = zeros(N, N)
    s = (s .== 0)*N + s
    d = (d .== 0)*N + d
for i in 1:n
    A[i,i] = 1.0 + A[s[i], d[i]]/2.0
        for j in (i+1):n
            if j > n break end
                A[i,j] = ( A[i, s[j]] + A[i,d[j]] )/2.0
                A[j,i] = A[i,j]
    end
    end
return(A[1:n, 1:n])
end

"""

	makePed(inputFile::String,userData)
Makes pedigree using PedigreeBase package

"""
function makePed(inputFile::String,userDataIDs)
	pedlist,idtable = read_ped(inputFile)
	
	perm,invp = find_ped_order(pedlist)
	permute_ped!(invp,pedlist,idtable)
	idtable = sort(idtable; byvalue=true)
	origIDs = [k for (k,v) in idtable if v in invp]

	issubset(userDataIDs,origIDs) || throw(ErrorException("ErrorException: phenotyed individuals are not a subset of pedigree"))

	f = get_inb(pedlist)
#       A = get_nrm(pedlist)
        Ainv = get_nrminv(pedlist, f)

	pedlist = DataFrame([origIDs invp pedlist'],:auto) #overwriting existing for memory
	rename!(pedlist,[:origID,:ID,:Sire,:Dam])
	return(pedlist,Ainv)
end


"""
	makeG(inputFile::String;method=1)
Makes genomic relationship matrix based on vanRaden method 1 (defult) or method 2
"""
function makeG(inputFile::String;method=1)
	thisM = CSV.read(inputFile,CSV.Tables.matrix,header=false)
	p = mean(thisM,dims=1)./2.0
	q = 1.0 .- p
        thisM .-= mean(thisM,dims=1) 
	if method==1 
		G = (thisM*thisM') ./ sum(2.0 .* p.*q)
	elseif method==2
		sqrt2pq = sqrt(2.0 .* p.*q)
		replace!(sqrt2pq,Inf=>0)
		thisM ./= sqrt2pq
		G = (thisM*thisM')./length(p)
	else error("enter a valid method")
	end
	G .+= Matrix(0.001*I,size(thisM,1),size(thisM,1))  
	return G	
end

"""
        makeG(inputData::Array{Float64,2};method=1)
Makes genomic relationship matrix based on vanRaden method 1 (defult) or method 2
"""
function makeG(thisM::Array{Float64,2};method=1)
        p = mean(thisM,dims=1)./2.0
        q = 1.0 .- p
        thisM .-= mean(thisM,dims=1)
        if method==1
                G = (thisM*thisM') ./ sum(2.0 .* p.*q)
        elseif method==2
                sqrt2pq = sqrt(2.0 .* p.*q)
                replace!(sqrt2pq,Inf=>0)
                thisM ./= sqrt2pq
                G = (thisM*thisM')./length(p)
        else error("enter a valid method")
        end
	G .+= Matrix(0.001*I,size(thisM,1),size(thisM,1))
        return G
end

#make regions
function prep2RegionData(outPutFolder,markerSet,mapFile,fixedRegSize)
    accRegion = 0
    accRegionVec = [0]
    SNPgroups = []
    mapData = CSV.read(mapFile,header=true,DataFrame)

    if fixedRegSize==99
        println("fixedRedSize $fixedRegSize")
        snpInfoFinal = mapData[!,[:snpID,:snpOrder,:chrID]]
	snpInfoFinal.groupID = snpInfoFinal.chrID
        accRegion    = length(unique(mapData[!,:chrID]))
        elseif fixedRegSize==9999
            snpInfoFinal = mapData[:,[:snpID,:snpOrder,:chrID]]
            snpInfoFinal.groupID  .= 1
            accRegion    = 1
        else
	
	snpInfoTemp = DataFrame(snpID=Vector{Any}(missing, 0),
                          groupID=Vector{Any}(missing, 0), copycols=false)
		
       	for c in unique(mapData[!,:chrID])
	    mapData = mapData[!,[:snpID,:snpOrder,:chrID]]	
            thisChr = mapData[mapData[!,:chrID] .== c,:]
            totLociChr = size(thisChr,1)
            TotRegions = ceil(Int,totLociChr/fixedRegSize)
            accRegion += TotRegions
            push!(accRegionVec, accRegion)
            tempGroups = sort(repeat(collect(accRegionVec[c]+1:accRegionVec[c+1]),fixedRegSize))

	tempInfo = DataFrame(snpID=Vector{Any}(missing, length(tempGroups)),
                          groupID=Vector{Any}(missing, length(tempGroups)), copycols=false)


            tempInfo[1:totLociChr,:snpID] = thisChr[!,:snpID]
            tempInfo[!,:groupID] = tempGroups
            dropmissing!(tempInfo)
            snpInfoTemp = vcat(snpInfoTemp,tempInfo)
            @printf("chr %.0f has %.0f groups \n", c, TotRegions)
	    grp = groupby(tempInfo, [:groupID])
	    println([size(x,1) for x in grp])
        end

	snpInfoFinal = deepcopy(mapData)
	snpInfoFinal.groupID = 	snpInfoTemp[!,:groupID]
	
        end  #ends if control flow
	CSV.write(outPutFolder*"/groupInfo_$(markerSet).txt",snpInfoFinal,delim='\t',header=true)
    for g in 1:accRegion
        push!(SNPgroups,searchsorted(snpInfoFinal[!,:groupID], g))
    end
    GC.gc()
    return SNPgroups
end

MatByMat = function(mat)
	mat'*mat
end

folderHandler = function(outFolder)
	if isdir(outFolder)==true
                println("Output folder $outFolder exists. Removing its content")
		rm(outFolder,force=true,recursive=true)
		sleep(1)
                run(`mkdir $outFolder`)

        else
                println("$outFolder has been created to store the MCMC output")
                run(`mkdir $outFolder`)
	end
end

myUnzip(d::Dict) = Dict(p.first => (;p.second...) for p in d)
myUnzip(d) = d



