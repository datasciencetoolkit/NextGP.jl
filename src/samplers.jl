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
			
				               		
        		#WRITE TO FILES
			if iter in these2Keep
				if isa(ySet,Symbol)
					inOut.outMCMC(outPut,"b_$ySet",hcat(vcat(b...)...))
					inOut.outMCMC(outPut,"varE_$ySet",hcat(varE[ySet]...))
				elseif isa(ySet,Tuple)
					inOut.outMCMC(outPut,"b_"*join(ySet, "_"),hcat(vcat(b...)...))
					inOut.outMCMC(outPut,"varE_"*join(ySet, "_"),hcat(varE[ySet]...))
				end

				[[inOut.outMCMC(outPut,"u$zSet$eSet",u[Z[zSet].pos]) for zSet in keys(Z) if (isa(zSet,Union{Symbol,Expr}) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single trait	
				[[inOut.outMCMC(outPut,"u"*join([String.(zSet)...])*"$eSet",hcat(u[Z[zSet].pos]...)) for zSet in keys(Z) if (isa(zSet,Tuple) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single-trait multiple comp
				
				#CHECK
				for zSet in keys(Z)
					inOut.outMCMC(outPut,"varU$zSet",hcat(reduce(hcat,varU[zSet])...))
				end
			
				[[inOut.outMCMC(outPut,"beta$mSet$eSet",beta[M[mSet].pos]) for mSet in keys(M) if (isa(mSet,Symbol) && in(mSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single trait	
				[[inOut.outMCMC(outPut,"delta$mSet$eSet",delta[M[mSet].pos]) for mSet in keys(M) if (isa(mSet,Symbol) && in(mSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single trait

				[[[inOut.outMCMC(outPut,"beta$m$eSet",beta[M[mSet].pos][[p],:]) for (p,m) in enumerate(mSet)] for mSet in keys(M) if (isa(mSet,Tuple{Vararg{Symbol}}) && in(mSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single-trait multiple comp
				#[[[inOut.outMCMC(outPut,"beta$m$eSet",beta[M[mSet].pos][[p],:]) for (p,m) in enumerate(mSet)] for mSet in keys(M) if (isa(mSet,Tuple{Vararg{Tuple{Vararg{Symbol}}}}) && in(mSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Tuple)] #multi-trait only one correlated comp	
				[[[inOut.outMCMC(outPut,"beta$m$(eSet[p])",beta[M[mSet].pos][[p],:]) for (p,m) in enumerate(mSet)] for mSet in keys(M) if (isa(mSet,Tuple{Vararg{Symbol}}))] for eSet in keys(E) if isa(eSet,Tuple)] #multi-trait only one correlated comp (:M1,M1) or (:M1,:M2...) but not ((:M1,:M2) and (:M3,:M4))
				[[[inOut.outMCMC(outPut,"delta$m$(eSet[p])",delta[M[mSet].pos][[p],:]) for (p,m) in enumerate(mSet)] for mSet in keys(M) if (isa(mSet,Tuple{Vararg{Symbol}}))] for eSet in keys(E) if isa(eSet,Tuple)] #multi-trait only one correlated comp (:M1,M1) or (:M1,:M2...) but not ((:M1,:M2) and (:M3,:M4))
				#[[[inOut.outMCMC(outPut,"delta$m$(eSet[p])",getindex.(M[mSet][:gammaHat],p)') for (p,m) in enumerate(mSet)] for mSet in keys(M) if (isa(mSet,Tuple{Vararg{Symbol}}))] for eSet in keys(E) if isa(eSet,Tuple)] #multi-trait only one correlated comp (:M1,M1) or (:M1,:M2...) but not ((:M1,:M2) and (:M3,:M4))

				#currently only one component with annotations
				[[inOut.outMCMC(outPut,"c$mSet$eSet",M[mSet][:c]) for mSet in keys(M) if (isa(mSet,Symbol) && in(mSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single trait	

				#for mSet in keys(M)
				#	if isa(mSet,Symbol)
				#		inOut.outMCMC(outPut,"beta$mSet",beta[M[mSet].pos])
				#		inOut.outMCMC(outPut,"delta$mSet",delta[M[mSet].pos])
				#		if in(M[mSet].method,["BayesB","BayesC","BayesR"])
				#			inOut.outMCMC(outPut,"pi$mSet",[M[mSet].piHat])	
				#		end
				#		if in(M[mSet].method,["BayesRCπ","BayesRCplus"])
				#			inOut.outMCMC(outPut,"pi$mSet",[vcat(M[mSet].piHat...)])
				#			inOut.outMCMC(outPut,"annot$mSet",[M[mSet].annotCat])	
				#		end
				#		if in(M[mSet].method,["BayesLV"])
				#			inOut.outMCMC(outPut,"c$mSet",[vcat(M[mSet].c...)])
				#			inOut.outMCMC(outPut,"varZeta$mSet",M[mSet].varZeta)
				#		end
				#	elseif isa(mSet,Tuple)
				#		for p in M[mSet].pos
				#			mSet2print = mSet[p]
				#			inOut.outMCMC(outPut,"beta$mSet2print",beta[p])	
				#		end
				#	end
                        	#end

				##scale and df is only one value/matrix array... adapt it here for speed
				for mSet in keys(M)
					if isa(mSet,Symbol)
						inOut.outMCMC(outPut,"var$(mSet)",hcat(reduce(hcat,varBeta[mSet])...))
						inOut.outMCMC(outPut,"scale$(mSet)",hcat(reduce(hcat,M[mSet].scale)...)) #df and scale is only one value array, later all should be converted from arrays to Float in all package
						inOut.outMCMC(outPut,"df$(mSet)",hcat(reduce(hcat,M[mSet].df)...))
					elseif isa(mSet,Tuple)
						inOut.outMCMC(outPut,"var_"*join(mSet, "_"),hcat(reduce(hcat,varBeta[mSet])...))
						inOut.outMCMC(outPut,"scale_"*join(mSet, "_"),hcat(reduce(hcat,M[mSet].scale)...))
						inOut.outMCMC(outPut,"df_"*join(mSet, "_"),hcat(reduce(hcat,M[mSet].df)...))
					end
				end
				
			end #end output freq.
		end #end ySet
	end #end iter
end #end function



#end #end module
