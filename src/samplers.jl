module samplers

using Distributions, LinearAlgebra
using StatsBase
using Printf

export runSampler

#main sampler
runSampler = function(Y,X,Z,varE,chainLength,burnIn,outputFreq) ##varE will be fixed for now
	
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #This is not really nFix, but the "blocks" of fixed effects
        nFix = length(X)

        #initial computations and settings
        
	#make iXpX
        iXpX = similar(X)
        for x in 1:nFix
                iXpX[x] = inv(X[x]'X[x])
        end
        #make b array
        b = Array{Array{Float64, 1},1}(undef,0)
        #counts columns per effect
        nColEachX = []
        for x in 1:nFix
                println(x)
                nCol = size(X[x],2)
                push!(b,fill(0.0,nCol))
                nColEachX = push!(nColEachX,nCol)
        end

        for iter in 1:chainLength
		#sample fixed effects
        	#always returns corrected Y and new b
        	sampleX!(X,b,iXpX,nFix,nColEachX,Y,varE)

        	#print
		if iter in these2Keep
	#		if onScreen==true
            			@printf("iter %s b %.2f \n", iter, vcat(b...)) #i always vectorize b. maybe better to make it vector initially
        #		end
		end
	end
end


#Sampling fixed effects
sampleX! = function(X,b,iXpX,nFix,nColEachX,ycorr,varE)
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

end
