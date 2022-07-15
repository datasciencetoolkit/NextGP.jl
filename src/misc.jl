using Printf
using DataFrames
using CSV
using PedigreeBase

"""
	makeA(s::Any, d::Any)
Makes pedigree-based relationship matrix.
adapted from http://morotalab.org/Mrode2005/relmat/createA.txt

"""

# adapted from http://morotalab.org/Mrode2005/relmat/createA.txt
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

	makePed!(inputFile::String,userData)
Makes pedigree using PedigreeBase package

"""

# using PedigreeBase package
function makePed(inputFile::String,userDataIDs)
	pedlist,idtable = read_ped(inputFile)
	
	issubset(userDataIDs,pedlist) || throw(ErrorException("ErrorException: phenotyed individuals are not a subset of pedigree")) 

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
                println("Output folder $outFolder exists. Removing it")
#                run(`rm -rf $outFolder/`)
		rm(outFolder,force=true,recursive=true)
        else
                println("$outFolder has been created to store the MCMC output")
                run(`mkdir $outFolder`)
	end
end






