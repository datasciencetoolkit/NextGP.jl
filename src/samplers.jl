module samplers

using Distributions, LinearAlgebra
using StatsBase
using Printf

include("outFiles.jl")

export runSampler

#main sampler
function runSampler(Y,X,Z,varE,chainLength,burnIn,outputFreq,Ai) ##varE will be fixed for now
	
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #This is not really nFix, but the "blocks" of fixed effects
	nData = length(Y)
        nFix  = length(X)
	nRand = length(Z)
	
        #initial computations and settings
        
	#make iXpX, Z', ZpZ
        iXpX = similar(X)
        for x in 1:nFix
                iXpX[x] = inv(X[x]'X[x])
        end
	Zp  = similar(Z')
	ZpZ = similar(Z)
	for z in 1:nRand
                ZpZ[z] = Z[z]'Z[z]
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
        for z in 1:nRand
                println(z)
                nCol = size(Z[z],2)
                push!(u,fill(0.0,nCol))
                nColEachZ = push!(nColEachZ,nCol)
        end


        for iter in 1:chainLength
		#sample fixed effects
        	#always returns corrected Y and new b
        	sampleX!(X,b,iXpX,nFix,nColEachX,Y,varE)
		sampleZ!(Ai,Zp,ZpZ,varE,varU,u,ycorr)
        	#print
		if iter in these2Keep
			IO.outMCMC(pwd(),vcat(b...)') ### currently no path is provided!!!!
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
		println("sampling from multi-variate normal")
		end
        	ycorr    .-= X[x]*b[x]
	end
end

#Sampling random effects
function sampleZ!(iMat,ZpMat,ZpZMat,varE,varU,u,ycorr)
	#block for each effect
	for z in 1:nRand
		λ = varE/varU	
	        ycorr .+= ZMat[z]*uVec[z]
	        Yi = ZpMat[z]*ycorr
		nCol = size(ZpZMat[z],2)
		uVec = deepcopy(u[z])
	        for i in 1:nCol
        	        uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
              		rhsU = Yi[i] - λ*dot(view(iMat,:,i),uVec)
                	lhsU = ZpZMat[i] + view(iMat,i,i)*λ
                	meanU = invLhsU*rhsU
                	uVec[i] = rand(Normal(meanU,sqrt(invLhsU*varE)))
        	end
		u[z] = uVec
        	ycorr .-= ZMat*uVec
	end
end


end
