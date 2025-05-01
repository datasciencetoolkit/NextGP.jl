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
include("types.jl")

include("functions.jl")
using .functions

export runSampler!

#main sampler
function runSampler!(ycorr,nData,E,X,b,Z,u,varU,M,beta,varBeta,delta,chainLength,burnIn,outputFreq,outPut)
		
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
		
		#sample residual variance
		for eSet in keys(E)
			varE[eSet] = sampleVarE!(eSet,E,varE,ycorr,nData)			
		end
		
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
#			println("running $(M[mSet].method) for $mSet")
			M[mSet].funct(mSet,M,beta,delta,ycorr,varE,varBeta)
		end
				               		
        	#print
		if iter in these2Keep
			inOut.outMCMC(outPut,"b",b') ### currently no path is provided!!!!
			inOut.outMCMC(outPut,"varE",varE)
			
			for zSet in keys(Z)
				if isa(zSet,Union{Expr,Symbol})
					inOut.outMCMC(outPut,"u$zSet",u[Z[zSet].pos])
				elseif isa(zSet,Tuple)
					pCounter = 1
					for p in zSet
						#zSet2print = zSet[p]
						inOut.outMCMC(outPut,"u$p",u[Z[zSet].pos][[pCounter],:])
						pCounter += 1
					end
				end
                        end

			for zSet in keys(Z)
				inOut.outMCMC(outPut,"varU$zSet",hcat(reduce(hcat,varU[zSet])...))
			end
			

			for mSet in keys(M)
				if isa(mSet,Symbol)
					inOut.outMCMC(outPut,"beta$mSet",beta[M[mSet].pos])
					inOut.outMCMC(outPut,"delta$mSet",delta[M[mSet].pos])
					if in(M[mSet].method,["BayesB","BayesC","BayesR"])
						inOut.outMCMC(outPut,"pi$mSet",[M[mSet].piHat])	
					end
					if in(M[mSet].method,["BayesRCπ","BayesRCplus"])
						inOut.outMCMC(outPut,"pi$mSet",[vcat(M[mSet].piHat...)])
						inOut.outMCMC(outPut,"annot$mSet",[M[mSet].annotCat])	
					end
					if in(M[mSet].method,["BayesLV"])
						inOut.outMCMC(outPut,"c$mSet",[vcat(M[mSet].c...)])
						inOut.outMCMC(outPut,"varZeta$mSet",M[mSet].varZeta)
					end
				elseif isa(mSet,Tuple)
					for p in M[mSet].pos
						mSet2print = mSet[p]
						inOut.outMCMC(outPut,"beta$mSet2print",beta[p])	
					end
				end
                        end

			for pSet in keys(M)
				inOut.outMCMC(outPut,"var$(pSet)",hcat(reduce(hcat,varBeta[pSet])...))
				inOut.outMCMC(outPut,"scale$(pSet)",hcat(reduce(hcat,M[pSet].scale)...))
			end
		end
	end
end



end
