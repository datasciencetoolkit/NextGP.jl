module samplers

using Distributions, LinearAlgebra
using StatsBase
using Printf
using CSV
using DataFrames

include("outFiles.jl")
include("misc.jl")

export runSampler

#main sampler
function runSampler(rowID,Y,X,Z,chainLength,burnIn,outputFreq,priorVCV,varM_prior,M,paths2maps,rS) ##varE will be fixed for now
	
	println("varM_prior $(varM_prior)")
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	nRand = length(Z)
	nData = length(Y)
	nMarkerSets = length(M)
	nMarkers    = [size(M[m],2) for m in 1:nMarkerSets]
	println("nMarkers: $nMarkers")

        #initial computations and settings
	ycorr = deepcopy(Y)
	        
	##make iXpX, Z', zpz (for uncor)
        iXpX = similar(X)
        for x in 1:nFix
                iXpX[x] = inv(X[x]'X[x])
        end

	Zp  = similar(Z) #similar(Z')
	zpz = Array{Array{Float64, 1},1}(undef,0)
	for z in 1:nRand
                push!(zpz,diag(Z[z]'Z[z]))
		Zp[z]  = Z[z]'
        end
	
        ##make b and u arrays
        b = Array{Array{Float64, 1},1}(undef,0)
        ##counts columns per effect
        nColEachX = []
        for x in 1:nFix
                println(x)
                nCol = size(X[x],2)
                push!(b,fill(0.0,nCol))
                nColEachX = push!(nColEachX,nCol)
        end

        u = Array{Array{Float64, 1},1}(undef,0)
        ##counts columns per effect
        nColEachZ = []
	##get priors per effect
	iVarStr = [] #inverses will be computed
	varU_prior = []
        for z in 1:nRand
                println(z)
                nCol = size(Z[z],2)
                push!(u,fill(0.0,nCol))
                nColEachZ = push!(nColEachZ,nCol)
		#var structures and priors
		if isempty(priorVCV[z][1])
			println("priorVCV $z is empty, an identity matrix will be used")
			push!(iVarStr,Matrix(1.0I,nCol,nCol))
		else 	push!(iVarStr,inv(priorVCV[z][1]))
		end
		push!(varU_prior,priorVCV[z][2])
        end
	println("prior variances $(varU_prior)")

	#set up for E	
	isempty(last(priorVCV,1)[1][1]) ? strE = Matrix(1.0I,nData,nData) : strE = last(priorVCV,1)[1][1]
	varE_prior = last(priorVCV,1)[1][2] #since last returns a tupple

	#parameters for priors
        dfE = 4.0
	dfDefault = 4.0
 
	dfM = Array{Float64}(undef,0)
        for m in 1:nMarkerSets
		size(varM_prior[m],1)>1 ? println("multivariate prior for marker set $m df=$(3+size(varM_prior[m],1))") : println("univariate prior for marker set $m df=$(3+size(varM_prior[m],1))")
                push!(dfM,3.0+size(varM_prior[m],1))
        end

	println("dfM $dfM")	
	       
	if varE_prior==0.0
		varE_prior  = 0.0005
       		scaleE     = 0.0005
        else
       		scaleE    = varE_prior*(dfE-2.0)/dfE    
   	end
	
	##no 0.0005 prior adapted here yet, also no correlated random effects
	scaleU = zeros(nRand)
	for z in 1:nRand
		scaleU[z] = varU_prior[z]*(dfDefault-2.0)/dfDefault
	end	

	scaleM = Array{Any,1}(undef,nMarkerSets)
	for m in 1:nMarkerSets
		nMComp = size(varM_prior[m],1)
                nMComp > 1 ? scaleM[m] = varM_prior[m].*(dfM[m]-nMComp-1.0)  : scaleM[m] = varM_prior[m]*(dfM[m]-2.0)/dfM[m] #I make float and array of float
        end


	#pre-computations using priors, not relevant for correlated random effects
   	νS_E = scaleE*dfE
	νS_U = zeros(nRand)
	for z in 1:nRand
                νS_U[z] = scaleU[z]*dfDefault
        end



	#ADD MARKERS
		# read map file and make regions
	regionArray =  Array{Array{UnitRange{Int64},1},1}(undef,0)
	for m in 1:nMarkerSets
		theseRegions = prep2RegionData(paths2maps[m],rS[m]) ###first data
		push!(regionArray,theseRegions)
	end
	println("size regionArray: $(length(regionArray))")
	println("size regionArray: $(length.(regionArray))")
	
	nRegions  = length.(regionArray) #per component

		#make mpm
#		Mp  = similar(M) #not needed coz I use BLAS.dot
       		mpm = Array{Array{Float64, 1},1}(undef,0)
       		for m in 1:nMarkerSets
               		push!(mpm,diag(M[m]'M[m])) #will not work for large matrices!!!!
#                	Mp[m]  = M[m]'
        	end

		#storage

	varU = varU_prior #for storage

	###FIXED FOR NOW
	mpmMat = 0

	beta = zeros(Float64,nMarkerSets,maximum(nMarkers)) #can allow unequal length! Remove tail zeros for printing....
#	vcovBeta = fill(Matrix(Diagonal(varM)),maximum(nRegions)) #can allow unequal length! Remove tail zeros for printing....

	varBeta = Array{Any,1}(undef,0)
        for m in 1:nMarkerSets
                push!(varBeta,fill(varM_prior[m],nRegions[m]))
        end

	#Start McMC
        for iter in 1:chainLength
		#sample residual variance
	       	varE = sampleVarE(νS_E,ycorr,dfE,nData)
		
		#sample fixed effects
        	#always returns corrected Y and new b
        	sampleX!(X,b,iXpX,nFix,nColEachX,ycorr,varE)
	
		#sample random effects
		# always returns corrected Y and new u
		sampleZ!(iVarStr,Z,Zp,zpz,nRand,varE,varU,u,ycorr)

		#sample variances
		sampleRanVar!(varU,nRand,νS_U,u,dfDefault,iVarStr)
		
		#sample marker effects
		sampleM!(M,beta,mpm,nMarkerSets,regionArray,ycorr,varE,varBeta)

		#sample marker variances
		

        	#print
		if iter in these2Keep
			IO.outMCMC(pwd(),"b",vcat(b...)') ### currently no path is provided!!!!
			IO.outMCMC(pwd(),"u",vcat(u...)')
			IO.outMCMC(pwd(),"varE",varE)
			IO.outMCMC(pwd(),"varU",varU')
	#		if onScreen==true
            			println("b, $(vcat(b...))") #i always vectorize b. maybe better to make it vector initially
        #		end
		end
	end
end


#Sampling fixed effects
function sampleX!(X,b,iXpX,nFix,nColEachX,ycorr,varE)
	#block for each effect 
	for x in 1:nFix
		ycorr    .+= X[x]*b[x]
        	rhs      = X[x]'*ycorr
                meanMu   = iXpX[x]*rhs
		if nColEachX[x] == 1
        		b[x] .= rand(Normal(meanMu[],sqrt((iXpX[x]*varE))[]))
		else b[x] .= rand(MvNormal(vec(meanMu),convert(Array,Symmetric(iXpX[x]*varE))))
		end
        	ycorr    .-= X[x]*b[x]
	end
end

#Sampling random effects
function sampleZ!(iStrMat,Zmat,ZpMat,zpzMat,nRand,varE,varU,u,ycorr)
	#block for each effect
	for z in 1:nRand
		uVec = deepcopy(u[z])
		iMat = iStrMat[z]
		tempzpz = zpzMat[z] ###added
		λz = varE/(varU[z])
	        ycorr .+= Zmat[z]*uVec		
	        Yi = ZpMat[z]*ycorr #computation of Z'ycorr for ALL  rhsU
		nCol = length(zpzMat[z])
	        for i in 1:nCol
        	        uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
			rhsU = Yi[i] - λz*dot(iMat[:,i],uVec)
                	lhsU = tempzpz[i] + (iMat[i,i]*λz)[1]
			invLhsU = 1.0/lhsU
                	meanU = invLhsU*rhsU
                	uVec[i] = rand(Normal(meanU,sqrt(invLhsU*varE)))
        	end
		u[z] = uVec
        	ycorr .-= Zmat[z]*uVec
	end
end


#Sampling marker effects
function sampleM!(MMat,beta,mpmMat,nMSet,regionsMat,ycorr,varE,varM)
        #for each marker set
        for mSet in 1:nMSet
		for r in 1:length(regionsMat[mSet]) #dont have to compute 1000000 times, take it out
			theseLoci = regionsMat[mSet][r]
			regionSize = length(theseLoci)
#			println("mSet: $mSet $regionSize $theseLoci")
			lambda = varE/(varM[mSet][r])
			println("mSet: $mSet reg size: $regionSize lambda: $lambda")
			for locus in theseLoci
				BLAS.axpy!(beta[mSet,locus],MMat[mSet][:,locus],ycorr)
				rhs = BLAS.dot(MMat[mSet][:,locus],ycorr)
				lhs = mpmMat[mSet][locus] + lambda
				meanBeta = lhs\rhs
				beta[mSet,locus] = sampleBeta(meanBeta, lhs, varE)
				BLAS.axpy!(-1.0*beta[mSet,locus],MMat[mSet][:,locus],ycorr)
			end
		end	
        end
end

#sample marker effects
function sampleBeta(meanBeta, lhs, varE)
    return rand(Normal(meanBeta,sqrt(lhs\varE)))
end

#sample random effects' variances
function sampleRanVar!(varU,nRand,νS_ranVar,effVec,df_ranVar,iStrMat)
	for z in 1:nRand
		n = size(iStrMat[z],2)
		varU[z] = (νS_ranVar[z] + effVec[z]'*iStrMat[z]*effVec[z])/rand(Chisq(df_ranVar + n))
	end
end

#Sample residual variance
function sampleVarE(νS_e,yCorVec,df_e,nRecords)
    return((νS_e + dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords)))
end



end
