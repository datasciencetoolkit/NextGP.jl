#module samplers


#using Distributions, LinearAlgebra
#using StatsBase
#using Printf
#using CSV
#using DataFrames
#using DataStructures
#using ProgressMeter
#using PrettyTables

#include("outFiles.jl")
#include("misc.jl")
#include("types.jl")

#include("functions.jl")
#using .functions


#export runSampler!

#main sampler
function runSampler!(modelInformation,ycorr,nData,E,varE,X,b,Z,u,varU,M,beta,varBeta,delta,chainLength,burnIn,outputFreq,outPut)
		
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
		
		#sample residual variance
		for (ySet,yModel) in modelInformation
			#sample error variance all at once!!!
			sampleVarE!(ySet,E,varE,ycorr,nData)
			
			#sample fixed effects
			isempty(X) ? nothing : [sampleX!(xSet,X,b,ycorr,varE,ySet) for xSet in keys(yModel) if isa(yModel[xSet],FixedEffect)] 

			#sample random effects
			#println("$ySet,$yModel")
			isempty(Z) ? nothing : [sampleZ!(zSet,Z,u,ycorr,varE,ySet,varU) for zSet in keys(yModel) if isa(yModel[zSet],RandomEffect)] 
	
			#sample marker effects and variances
			isempty(M) ? nothing : [M[mSet].funct(mSet,M,beta,delta,ycorr,varE,varBeta,ySet) for mSet in keys(yModel) if isa(yModel[mSet],RandomEffect)] 
			
#			for mSet in keys(M)
#				M[mSet].funct(mSet,M,beta,delta,ycorr,varE,varBeta)
#			end
				               		
        		#WRITE TO FILES
			if iter in these2Keep
				inOut.outMCMC(outPut,"b_$ySet",hcat(vcat(b...)...))
				inOut.outMCMC(outPut,"varE_$ySet",hcat(varE[ySet]...))

				[[inOut.outMCMC(outPut,"u$zSet$eSet",u[Z[zSet].pos]) for zSet in keys(Z) if (isa(zSet,Union{Symbol,Expr}) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single trait	
				[[inOut.outMCMC(outPut,"u"*join([String.(zSet)...])*"$eSet",hcat(u[Z[zSet].pos]...)) for zSet in keys(Z) if (isa(zSet,Tuple) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single-trait multiple comp
				
				#CHECK
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



#end #end module
