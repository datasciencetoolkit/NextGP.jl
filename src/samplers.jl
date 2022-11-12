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
function runSampler!(ycorr,nData,dfE,scaleE,X,b,Z,iVarStr,Zp,zpz,uKeyPos,uKeyPos4Print,nColEachZ,u,varU,scaleZ,dfZ,M,beta,varBeta,BayesX,rhsZ,chainLength,burnIn,outputFreq,outPut)
		
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
	
		#sample residual variance
	       	varE = sampleVarE(dfE,scaleE,ycorr,nData)
		
		#sample fixed effects

		for xSet in keys(X)
			sampleX!(X[xSet],b,ycorr,varE)
		end
	
		#sample random effects
	        sampleZandZVar!(iVarStr,Z,Zp,u,zpz,uKeyPos,nColEachZ,ycorr,varE,varU,scaleZ,dfZ)	

		#sample marker effects and variances
	
		
		for mSet in keys(M)
			BayesX[mSet](mSet,M,beta,ycorr,varE,varBeta)
		end
				               		
        	#print
		if iter in these2Keep
			IO.outMCMC(outPut,"b",b') ### currently no path is provided!!!!
			IO.outMCMC(outPut,"varE",varE)
			
			for zSet in keys(uKeyPos4Print)
                                IO.outMCMC(outPut,"u$(uKeyPos4Print[zSet])",u[uKeyPos4Print[zSet],1:nColEachZ[zSet]]')
                        end
			for pSet in keys(zpz)
				IO.outMCMC(outPut,"varU$(uKeyPos[pSet])",varU[pSet]) #join values for multivariate in uKeyPos[pSet])
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
