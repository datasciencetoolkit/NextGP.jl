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

include("functions.jl")
using .functions

include("prepMatVec.jl")

export runSampler!

#main sampler
function runSampler!(modelInformation,ycorr,nData,E,varE,X,b,Z,u,varU,M,beta,varBeta,delta,chainLength,burnIn,outputFreq,outPut)
		
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
		
		#sample residual variance
		for (ySet,yModel) in modelInformation

			####NEED TO ADD VARE HERE
			
			#sample fixed effects

			#for xSet in yModel #keys(X)
			#	println("sampling now $xSet")
			#	sampleX!(xSet,X,b,ycorr,varE,ySet)
			#end

			#[sampleX!(xSet,X,b,ycorr,varE,ySet) for xSet in keys(yModel) if isa(yModel[xSet],DataTerm)]
			#[println("sampling $xSet") for xSet in keys(yModel) if isa(yModel[xSet],DataTerm)]
			gh = [typeof(yModel[xSet]) for xSet in keys(yModel) if typeof(yModel[xSet]) <: prepMatVec.FixedEffect]
			println("GH: $gh")
			
			#sample random effects and variances
			for zSet in keys(Z)
	        		sampleZ!(zSet,Z,u,ycorr,varE,varU)	
			end
		
	
			#sample marker effects and variances
			for mSet in keys(M)
#				println("running $(M[mSet].method) for $mSet")
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
			end #end output freq.
		end #end ySet
	end #end iter
end #end function



end #end module
