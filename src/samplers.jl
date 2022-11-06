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
function runSampler!(ycorr,nData,dfE,scaleE,X,iXpX,XKeyPos,b,Z,iVarStr,Zp,zpz,uKeyPos,uKeyPos4Print,nColEachZ,u,varU,scaleZ,dfZ,M,Mp,mpm,BetaKeyPos,BetaKeyPos4Print,beta,regionArray,nRegions,varBeta,scaleM,dfM,BayesX,chainLength,burnIn,outputFreq,outPut)
		
	println("BayesX in sampler: $BayesX")
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
	
		#sample residual variance
	       	varE = sampleVarE(dfE,scaleE,ycorr,nData)
		
		#sample fixed effects

		for xSet in keys(iXpX)
			sampleX!(xSet,X[xSet],b,iXpX[xSet],XKeyPos[xSet],ycorr,varE)
		end
	
		#sample random effects
	        sampleZandZVar!(iVarStr,Z,Zp,u,zpz,uKeyPos,nColEachZ,ycorr,varE,varU,scaleZ,dfZ)	

		#sample marker effects and variances
	
		
		for mSet in keys(mpm)
#			sampleBayesPR!(mSet,M[mSet],Mp[mSet],beta,mpm[mSet],BetaKeyPos[mSet],regionArray[mSet],nRegions[mSet],ycorr,varE,varBeta,scaleM[mSet],dfM[mSet])
			println("running: $(BayesX[mSet])")
			sampleBayesPR!(mSet,M[mSet],Mp[mSet],beta,mpm[mSet],BetaKeyPos[mSet],regionArray[mSet],nRegions[mSet],ycorr,varE,varBeta,scaleM[mSet],dfM[mSet])
			method = BayesX[mSet]
			println("running: $(BayesX[mSet])")
			functions.method(mSet,M[mSet],Mp[mSet],beta,mpm[mSet],BetaKeyPos[mSet],regionArray[mSet],nRegions[mSet],ycorr,varE,varBeta,scaleM[mSet],dfM[mSet])
		end
               		
        	#print
		if iter in these2Keep
			IO.outMCMC(outPut,"b",vcat(b...)') ### currently no path is provided!!!!
			IO.outMCMC(outPut,"varE",varE)
			
			for zSet in keys(uKeyPos4Print)
                                IO.outMCMC(outPut,"u$(uKeyPos4Print[zSet])",u[uKeyPos4Print[zSet],1:nColEachZ[zSet]]')
                        end
			for pSet in keys(zpz)
				IO.outMCMC(outPut,"varU$(uKeyPos[pSet])",varU[pSet]) #join values for multivariate in uKeyPos[pSet])
			end
			

			for mSet in keys(BetaKeyPos4Print)
				IO.outMCMC(outPut,"beta$mSet",beta[BetaKeyPos4Print[mSet]])
                        end

			for pSet in keys(mpm)
				IO.outMCMC(outPut,"var$(pSet)",hcat(reduce(hcat,varBeta[pSet])...))
			end
		end
	end
end

end
