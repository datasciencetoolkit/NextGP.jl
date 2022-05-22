module samplers

using Distributions, LinearAlgebra
using StatsBase
using Printf

include("outFiles.jl")

export runSampler

#main sampler
function runSampler(rowID,Y,X,Z,chainLength,burnIn,outputFreq,priorVCV) ##varE will be fixed for now
	

	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	nRand = length(Z)
	nData = length(Y)

        #initial computations and settings
	ycorr = deepcopy(Y)
	        
	#make iXpX, Z', zpz (for uncor)
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
	
        #make b and u arrays
        b = Array{Array{Float64, 1},1}(undef,0)
        #counts columns per effect
        nColEachX = []
        for x in 1:nFix
                println(x)
                nCol = size(X[x],2)
                push!(b,fill(0.0,nCol))
                nColEachX = push!(nColEachX,nCol)
        end

        u = Array{Array{Float64, 1},1}(undef,0)
        #counts columns per effect
        nColEachZ = []
	#get priors per effect
	iVarStr = [] #inverses will be computed
	varU = []
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
			println("USED MAT: $(priorVCV[z][1])")
		end
		push!(varU,priorVCV[z][2])
        end
	println("prior variances $(varU)")

	#set up for E	
	#variances are gonna be priors, but fixed now!
        ##############################################
	isempty(last(priorVCV,1)[1][1]) ? strE = Matrix(1.0I,nData,nData) : strE = last(priorVCV,1)[1][1]
#	strE = Matrix(1.0I,nData,nData)
	varE = last(priorVCV,1)[1][2] #since last returns a tupple
	println("structure for E: $strE")
	println("prior for E: $varE")

        for iter in 1:chainLength
		#sample fixed effects
        	#always returns corrected Y and new b
        	sampleX!(X,b,iXpX,nFix,nColEachX,ycorr,varE)
	
		#sample random effects
		# always returns corrected Y and new u
		sampleZ!(iVarStr,Z,Zp,zpz,nRand,varE,varU,u,ycorr)
        	#print
		if iter in these2Keep
			IO.outMCMC(pwd(),"b",vcat(b...)') ### currently no path is provided!!!!
			IO.outMCMC(pwd(),"u",vcat(u...)')
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
			println("sampling from uni-variate normal")
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


end
