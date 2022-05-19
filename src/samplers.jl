module samplers

using Distributions, LinearAlgebra
using StatsBase

export runSampler

#main sampler
runSampler = function(Y,X,Z,varE,nChain) ##varE will be fixed for now
        nFix = length(X)

        #initial computations and settings
        
	#make iXpX
        iXpX = similar(X)
        for x in 1:nFix
                iXpX[x] = inv(X[x]'X[x] + Matrix(0.001I,sum(size.(X,2)),sum(size.(X,2))))
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

        for i in 1:nChain
		#sample fixed effects
        	#always returns corrected Y and new b
        	sampleX!(X,b,iXpX,nFix,nColEachX,Y,varE)

        	#print
        	println("sampled b: $(vcat(b...))")
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
