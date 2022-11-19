module samplers


using Distributions, LinearAlgebra
using StatsBase
using Printf
using CSV
using DataFrames
using DataStructures
using ProgressMeter
using PrettyTables

include("outFiles.jl")
include("misc.jl")
include("runTime.jl")

include("functions.jl")
using .functions

export runSampler!

#main sampler
function runSampler!(ycorr,nData,dfE,scaleE,X,b,Z,u,varU,M,beta,varBeta,BayesX,chainLength,burnIn,outputFreq,outPut)
		
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
	
		#sample residual variance
	       	varE = sampleVarE(dfE,scaleE,ycorr,nData)
		
		#sample fixed effects

		for xSet in keys(X)
			sampleX!(xSet,X,b,ycorr,varE)
		end
	
		#sample random effects and variances
		for zSet in keys(Z)
	        	functions.sampleZ!(zSet,Z,u,ycorr,varE,varU)	
		end
		
	
		#sample marker effects and variances
		for mSet in keys(M)
			BayesX[mSet](mSet,M,beta,ycorr,varE,varBeta)
		end
				               		
        	#print
		if iter in these2Keep
			IO.outMCMC(outPut,"b",b') ### currently no path is provided!!!!
			IO.outMCMC(outPut,"varE",varE)
			
			for zSet in keys(Z)
				if isa(zSet,Union{Expr,Symbol})
					IO.outMCMC(outPut,"u$zSet",u[Z[zSet].pos])
				elseif isa(zSet,Tuple)
					pCounter = 1
					for p in zSet
						#zSet2print = zSet[p]
						println("printing $p to $(Z[zSet].pos) and $pCounter")
						IO.outMCMC(outPut,"u$p",u[Z[zSet].pos][[pCounter],:])
						pCounter += 1
					end
				end
                        end

			for zSet in keys(Z)
				IO.outMCMC(outPut,"varU$zSet",hcat(reduce(hcat,varU[zSet])...))
			end
			

			for mSet in keys(M)
				if isa(mSet,Symbol)
					IO.outMCMC(outPut,"beta$mSet",beta[M[mSet].pos])
				elseif isa(mSet,Tuple)
					for p in M[mSet].pos
						mSet2print = mSet[p]
						IO.outMCMC(outPut,"beta$mSet2print",beta[p])	
					end
				end
                        end

			for pSet in keys(M)
				IO.outMCMC(outPut,"var$(pSet)",hcat(reduce(hcat,varBeta[pSet])...))
			end
		end
	end
end

end
