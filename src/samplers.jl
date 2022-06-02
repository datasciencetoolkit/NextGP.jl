module samplers

using Distributions, LinearAlgebra
using StatsBase
using Printf

include("outFiles.jl")

export runSampler

#main sampler
function runSampler(rowID,Y,X,Z,chainLength,burnIn,outputFreq,priorVCV,M,map,rS) ##varE will be fixed for now
	

	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	nRand = length(Z)
	nData = length(Y)
	nMarkerSets = length(M)

        #initial computations and settings
	ycorr = deepcopy(Y)
	        
	##make iXpX, Z', zpz (for uncor)
        iXpX = similar(X)
        for x in 1:nFix
                iXpX[x] = inv(X[x]'X[x])
        end

	Zp  = similar(Z')
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
        dfE = 4
	dfDefault = 4 
	       
	if varE_prior==0.0
		varE_prior  = 0.0005
       		scaleE     = 0.0005
        else
       		scaleE    = varE_prior*(dfE-2.0)/dfE    
   	end
	
	##no 0.0005 prior adapted here yet
	scaleU = zeros(nRand)
	for z in 1:nRand
		scaleU[z] = varU_prior[z]*(dfDefault-2.0)/dfDefault
	end	


	#pre-computations using priors
   	νS_E = scaleE*dfE
	νS_U = zeros(nRand)
	for z in 1:nRand
                νS_U[z] = scaleU[z]*dfDefault
        end 


	#FIXED RANDOM EFFECTS' VARIANCES
	varU = varU_prior #for storage


	#ADD MARKERS
		# read map file and make regions
		
	

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
