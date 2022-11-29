module mme


using Distributions, LinearAlgebra
using StatsBase
using Printf
using CSV
using DataFrames
using DataStructures
using PrettyTables

include("outFiles.jl")
include("misc.jl")
include("runTime.jl")

include("functions.jl")
using .functions

export getMME!

#define type for priorVCV to include Expression :(1|Dam)  or Symbol (:Dam)
ExprOrSymbol = Union{Expr,Symbol,Tuple}

#set up (co)variance structures
function setVarCovStr!(zSet::ExprOrSymbol,Z::Dict,priorVCV,varU_prior::Dict)
	if haskey(priorVCV,zSet)	
		if ismissing(priorVCV[zSet].str) || priorVCV[zSet].str=="I" 
			printstyled("prior var-cov structure for $zSet is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
			Z[zSet][:iVarStr] = Matrix(1.0I,Z[zSet][:dims][2],Z[zSet][:dims][2])
		elseif priorVCV[zSet].str=="A"
			printstyled("prior var-cov structure for $zSet is A. Computed A matrix (from pedigree file) will be used\n"; color = :green)
			isa(zSet,Tuple) ? Z[zSet][:iVarStr] = Z[zSet[1]][:iVarStr] : Z[zSet][:iVarStr] = Z[zSet][:iVarStr]
		elseif priorVCV[zSet].str=="G"
                        printstyled("prior var-cov structure for $zSet is G. Computed G matrix will be used\n"; color = :green)
			isa(zSet,Tuple) ? Z[zSet][:iVarStr] = Z[zSet[1]][:iVarStr] : Z[zSet][:iVarStr] = Z[zSet][:iVarStr]
		else 	Z[zSet][:iVarStr] = inv(priorVCV[zSet].str)
		end
		varU_prior[zSet] = priorVCV[zSet].v
	else	
		printstyled("prior var-cov for $zSet is empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		varU_prior[zSet] = 100
		priorVCV[zSet] = Random("I",0,100)
		Z[zSet][:iVarStr] = Matrix(1.0I,Z[zSet][:dims][2],Z[zSet][:dims][2])
	end
end


#main sampler
function getMME!(Y,X,Z,M,blocks,priorVCV,summaryStat,outPut)
		
        #some info
	nRand = length(Z)
	nData = length(Y)
	nMarkerSets = length(M)
        #initial computations and settings
	ycorr = deepcopy(Y)

	priorVCV = convert(Dict{ExprOrSymbol, Any},priorVCV)	
	
	### X and b	
	
	for blk in blocks
		X[blk] = Dict{Symbol, Any}()
		X[blk][:data] = hcat(getindex.(getindex.(Ref(X), blk),:data)...)
		X[blk][:levels] = hcat(vcat(getindex.(getindex.(Ref(X), blk),:levels)...)...)
		X[blk][:nCol] = sum(getindex.(getindex.(Ref(X), blk),:nCol))
		X[blk][:method] = first(getindex.(getindex.(Ref(X), blk),:method))
		for d in blk
			delete!(X,d)
		end
	end

	#==BLOCK FIXED EFFECTS.
	Order of blocks is as definde by the user
	Order of variables within blocks is always the same as in the model definition, not defined by the user in each block.
	==#
	
	
	#Positions of parameters for each variable and blocks for speed. b is a column vector.
	countXCol = 0
	for xSet in keys(X)
		nCol = X[xSet][:nCol]
                X[xSet][:pos] = (countXCol+1):(countXCol+nCol)
		countXCol += nCol
        end
	countXCol = 0


	b = zeros(sum(getindex.(getindex.(Ref(X), keys(X)),:nCol)))


	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	
        for xSet in keys(X)
#		ixpx inverse taken later
		X[xSet][:ixpx] = X[xSet][:data]'X[xSet][:data]
		X[xSet][:rhs] = zeros(X[xSet][:nCol])

                if xSet in keys(summaryStat)
	  		summaryStat[xSet].v == Array{Float64,1} ? X[xSet][:ixpx] .+= inv.(summaryStat[xSet].v) : X[xSet][:ixpx] .+= inv.(diag(summaryStat[xSet].v))
			summaryStat[xSet].v == Array{Float64,1} ? X[xSet][:rhs] .= inv.(summaryStat[xSet].v) .* (summaryStat[xSet].m)  : X[xSet][:rhs] .= inv.(diag(summaryStat[xSet].v)) .* (summaryStat[xSet].m)
                end

		if isa(X[xSet][:ixpx],Matrix{Float64}) 
			X[xSet][:ixpx] += Matrix(I*minimum(abs.(diag(X[xSet][:ixpx])./10000)),size(X[xSet][:ixpx]))
		end
               	X[xSet][:ixpx] = inv(X[xSet][:ixpx])
        end
	        
      

	#set up for E.
						
	#no inverse implemented yet!
	if haskey(priorVCV,:e)	
		if ismissing(priorVCV[:e].str) || priorVCV[:e].str=="I" 
				printstyled("prior var-cov structure for \"e\" is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				strE = Matrix(1.0I,nData,nData)
				priorVCV[:e] = Random("I",priorVCV[:e].m,priorVCV[:e].v)
		elseif priorVCV[:e].str=="D"
				strE = D ##no inverse  yet
				error("var-cov structure \"D\" has not been implemented yet")
				printstyled("prior var-cov structure for \"e\" is \"D\". User provided \"D\" matrix (d_ii = 1/w_ii) will be used\n"; color = :green)
		else 
				error("provide a valid prior var-cov structure (\"I\", \"D\" or leave it empty \"[]\") for \"e\" ")
		end
	else	
		printstyled("prior var-cov for \"e\" is fully  empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		strE = Matrix(1.0I,nData,nData)
		#just add to priors
		priorVCV[:e] = Random("I",0,100)
	end
								
	#parameters for priors
        dfE = 4.0
 	       
	if priorVCV[:e].v==0.0
		priorVCV[:e].v  = 0.0005
       		scaleE     = 0.0005
        else
       		scaleE    = priorVCV[:e].v*(dfE-2.0)/dfE    
   	end


	#### New u
	
	#key positions for each effect in u, for speed. Order of matrices in Z are preserved here.
					
#        for zSet in keys(Z)
#		pos = findall(x->x==zSet, collect(keys(Z)))[]
#                Z[zSet][:pos] = pos
#		println("$zSet : $pos")
#        end
	
	u = []
	varU_prior = Dict{Any,Any}()

	#matrices are ready
	
	posZcounter = 0
	for zSet in keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV)
		#symbol :ID or expression :(1|ID)
		if (isa(zSet,Symbol) || isa(zSet,Expr)) && in(zSet,keys(Z))
			posZcounter += 1
			Z[zSet][:pos] = posZcounter
			tempzpz = []
			nowZ = Z[zSet][:data]
			setVarCovStr!(zSet,Z,priorVCV,varU_prior)
			for c in eachcol(nowZ)
				push!(tempzpz,c'c)					
				# push!(tempzpz,BLAS.dot(c,c))
			end
			Z[zSet][:zpz] = tempzpz
			Z[zSet][:rhs] = zeros(size(Z[zSet][:data],2))
                        if zSet in keys(summaryStat)
                                summaryStat[zSet].v == Array{Float64,1} ? zpz[zSet] .+= inv.(summaryStat[zSet].v) : zpz[zSet] .+= inv.(diag(summaryStat[zSet].v))
                                summaryStat[zSet].v == Array{Float64,1} ? rhsZ[zSet] .= inv.(summaryStat[zSet].v) .* (summaryStat[zSet].m)  : rhsZ[zSet] .= inv.(diag(summaryStat[zSet].v)) .* (summaryStat[zSet].m)
                        end
			Z[zSet][:Zp]  = transpose(Z[zSet][:data])
			u = push!(u,zeros(Float64,1,size(Z[zSet][:data],2)))
		#tuple of symbols (:ID,:Dam)
		elseif (isa(zSet,Tuple{Vararg{Symbol}})) && all((in).(zSet,Ref(keys(Z)))) #if all elements are available # all([zSet .in Ref(keys(Z))])
			Z[zSet] = Dict{Symbol, Any}()
			Z[pSet][:pos] = collect((posZcounter+1):(posZcounter+length(pSet)))
			posZcounter += length(pSet)
			#posZcounter += 1
			#Z[zSet][:pos] = posZcounter
			u = push!(u,zeros(Float64,length(zSet),size(Z[zSet[1]][:data],2)))
			Z[zSet][:levels] = first(getindex.(getindex.(Ref(Z),zSet),:levels))
			println("Z: $Z")
			tempZ = hcat.(eachcol.(getindex.(getindex.(Ref(Z), zSet),:data))...)
			#same Z for all components in a single-trait model get only first column! Z[zSet][:data] = getindex.(tempZ,:,1)
			Z[zSet][:data] = tempZ
			setVarCovStr!(zSet,Z,priorVCV,varU_prior)
			Z[zSet][:str] = Z[zSet[1]][:str] 
			for d in zSet
                       		delete!(Z,d)
               		end
			Z[zSet][:zpz]  = MatByMat.(tempZ)
			Z[zSet][:Zp]   = transpose.(tempZ)
			tempZ = 0
			#lhs is already zero as only mpm + "nothing" is  given
			#rhs is for now only for convenience
			Z[zSet][:rhs] = [zeros(length(zSet)) for i in 1:length(Z[zSet][:levels])]
			if zSet in keys(summaryStat)
				error("Not available to use summary statistics in correlated effects")
                        	#SummaryStat[pSet].v == Array{Float64,1} ? zpz[pSet] += inv.(SummaryStat[pSet].v) : zpz[pSet] += inv.(diag(SummaryStat[pSet].v))
                        	#SummaryStat[pSet].v == Array{Float64,1} ? rhsZ[pSet] = inv.(SummaryStat[pSet].v) .* (SummaryStat[pSet].m)  : rhsZ[pSet] = inv.(diag(SummaryStat[pSet].v)) .* (SummaryStat[pSet].m)
                	end
#			Z[zSet][:Zp]  = transpose.(tempZ)
			tempZ = 0
		end
	end
	
	for zSet in collect(keys(Z))[(!in).(keys(Z),Ref(keys(priorVCV)))]
		posZcounter += 1
		Z[zSet][:pos] = posZcounter
		printstyled("No prior was provided for $pSet, but it was not included in the data. It will be made uncorrelated with default priors\n"; color = :green)		
		tempzpz = []
		nowZ = Z[zSet][:data]
		for c in eachcol(nowZ)
			push!(tempzpz,c'c)					
		end
		Z[zSet][:Zp]  = transpose(Z[zSet])						
		Z[zSet][:zpz] = tempzpz
		Z[zSet][:rhs] = zeros(size(Z[zSet][:dims],2))
		if zSet in keys(summaryStat)
                	summaryStat[zSet].v == Array{Float64,1} ? zpz[zSet] .+= inv.(summaryStat[zSet].v) : zpz[zSet] .+= inv.(diag(summaryStat[zSet].v))
                        summaryStat[zSet].v == Array{Float64,1} ? rhsZ[zSet] .= inv.(summaryStat[zSet].v) .* (summaryStat[zSet].m)  : rhsZ[zSet] .= inv.(diag(summaryStat[zSet].v)) .* (summaryStat[zSet].m)
                end
		setVarCovStr!(zSet,Z,priorVCV,varU_prior)
	end
																		
	#df, shape, scale...															
	
	for zSet ∈ keys(Z)
		Z[zSet][:df] = 3.0+size(priorVCV[zSet].v,1)
	end
																
        for zSet in keys(Z)
                nZComp = size(priorVCV[zSet].v,1)
		#priorVCV[zSet].v is a temporary solution
		nZComp > 1 ? Z[zSet][:scale] = priorVCV[zSet].v .* (Z[zSet][:df]-nZComp-1.0)  : Z[zSet][:scale] = priorVCV[zSet].v * (Z[zSet][:df]-2.0)/Z[zSet][:df] #I make float and array of float														
        end


												
        ####
																					

	#ADD MARKERS
	# read map file and make regions
																		
	############priorVCV cannot be empty for markers, currently!!																	

	beta = []

	#make mpm
	posMcounter = 0
	for pSet in keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV)
		#symbol :M1 or expression
		if isa(pSet,Symbol) && in(pSet,keys(M))
			posMcounter += 1
			M[pSet][:pos] = posMcounter
			tempmpm = []
			nowM = M[pSet][:data]
			for c in eachcol(nowM)
				push!(tempmpm,BLAS.dot(c,c))
			end
			M[pSet][:mpm] = tempmpm
			M[pSet][:rhs] = zeros(M[pSet][:dims][2])
			if pSet in keys(summaryStat)
                       		summaryStat[pSet].v == Array{Float64,1} ? M[pSet][:mpm] .+= inv.(summaryStat[pSet].v) : M[pSet][:mpm] .+= inv.(diag(summaryStat[pSet].v))
				summaryStat[pSet].v == Array{Float64,1} ? M[pSet][:rhs] .= inv.(summaryStat[pSet].v) .* (summaryStat[pSet].m)  : M[pSet][:rhs] .= inv.(diag(summaryStat[pSet].v)) .* (summaryStat[pSet].m)
			end
			M[pSet][:Mp] = []
			theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],priorVCV[pSet].r)
		        M[pSet][:regionArray] = theseRegions
			M[pSet][:nRegions] = length(theseRegions)
			beta = push!(beta,zeros(Float64,1,M[pSet][:dims][2]))	
		#tuple of symbols (:M1,:M2)
		elseif (isa(pSet,Tuple{Vararg{Symbol}})) && all((in).(pSet,Ref(keys(M)))) #if all elements are available # all([pSet .in Ref(keys(M))])
			M[pSet] = Dict{Symbol, Any}()
			M[pSet][:pos] = collect((posMcounter+1):(posMcounter+length(pSet)))
			posMcounter += length(pSet)
			maps = getindex.(getindex.(Ref(M),pSet),:map)
			(length(maps)==0 || all( ==(maps[1]), maps)) == true ? M[pSet][:map] = first(maps) : error("correlated marker sets must have the same map file!")
#			M[pSet][:pos] = vcat(getindex.(getindex.(Ref(M), pSet),:pos)...)
			M[pSet][:levels] = first(getindex.(getindex.(Ref(M),pSet),:levels))
			tempM = hcat.(eachcol.(getindex.(getindex.(Ref(M), pSet),:data))...)
			M[pSet][:data] = tempM
			for d in pSet
				beta = push!(beta,zeros(Float64,1,M[d][:dims][2]))
                       		delete!(M,d)
               		end
			M[pSet][:mpm] = MatByMat.(tempM)
			M[pSet][:Mp]  = transpose.(tempM)
			tempM = 0
			#lhs is already zero as only mpm + "nothing" is  given
			#rhs is for now only for convenience
			M[pSet][:rhs] = [zeros(length(pSet)) for i in 1:length(M[pSet][:levels])]
			if pSet in keys(summaryStat)
				error("Not available to use summary statistics in correlated effects")
                                #SummaryStat[pSet].v == Array{Float64,1} ? mpm[pSet] += (1.0 ./ SummaryStat[pSet].v) : mpm[pSet] += inv.(diag(SummaryStat[pSet].v))
 	                end
			theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],priorVCV[pSet].r)
			M[pSet][:regionArray] = theseRegions
			M[pSet][:nRegions] = length(theseRegions)
		end
	end
	
	for pSet in collect(keys(M))[(!in).(keys(M),Ref(keys(priorVCV)))]
		posMcounter += 1
		M[pSet][:pos] = posMcounter
		printstyled("No prior was provided for $pSet, but it was included in the data. It will be made uncorrelated with default priors and region size 9999 (WG)\n"; color = :green)		
		tempmpm = []
		nowM = M[pSet][:data]
		for c in eachcol(nowM)
			push!(tempmpm,BLAS.dot(c,c))
		end
		M[pSet][:mpm] = tempmpm
		M[pSet][:rhs] = zeros(M[pSet][:dims][2])
                if pSet in keys(summaryStat)
			summaryStat[pSet].v == Array{Float64,1} ? M[pSet][:mpm] .+= inv.(summaryStat[pSet].v) : M[pSet][:mpm] .+= inv.(diag(summaryStat[pSet].v))
                        summaryStat[pSet].v == Array{Float64,1} ? M[pSet][:rhs] .= inv.(summaryStat[pSet].v) .* (summaryStat[pSet].m)  : M[pSet][:rhs] .= inv.(diag(summaryStat[pSet].v)) .* (summaryStat[pSet].m)
                end
		theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],9999)
		M[pSet][:regionArray] = theseRegions
		M[pSet][:nRegions] = length(theseRegions)
	end

	
	for mSet ∈ keys(M)
		M[mSet][:df] = 3.0+size(priorVCV[mSet].v,1)
	end


        for mSet in keys(M)
                nMComp = size(priorVCV[mSet].v,1)
                nMComp > 1 ? M[mSet][:scale] = priorVCV[mSet].v .* (M[mSet][:df]-nMComp-1.0)  : M[mSet][:scale] = priorVCV[mSet].v * (M[mSet][:df]-2.0)/(M[mSet][:df]) #I make float and array of float
        end
	
	
	#storage

	varU = deepcopy(varU_prior) #for storage

	varBeta = Dict{Union{Symbol,Tuple{Vararg{Symbol}}},Any}()
        for mSet in keys(M)
                varBeta[mSet] = [priorVCV[mSet].v for i in 1:M[mSet][:nRegions]] #later, direct reference to key when varM_prior is a dictionary
        end

	#summarize analysis
	summarize = DataFrame(Effect=Any[],Type=Any[],Str=Any[],df=Any[],scale=Any[])
	
	for zSet in keys(Z)
		push!(summarize,[zSet,"Random",Z[zSet][:str],Z[zSet][:df],Z[zSet][:scale]])
	end


	###Bayesian Alphabet methods
	BayesX = Dict{Union{Symbol,Tuple{Vararg{Symbol}}},Any}()

	for mSet in keys(M)
		if mSet ∈ keys(priorVCV)
			priorVCV[mSet].name == "BayesPR" ? BayesX[mSet] = sampleBayesPR! : nothing
			#BayesX[mSet] = typeof(priorVCV[mSet])
			str = "$(M[mSet][:nRegions]) block(s)"
			#value = priorVCV[mSet].v
		else #### later, handel this above, when dealing with priorVCV is allowed to be empty
			BayesX[mSet] = BayesPR #with region size 9999
			str = "WG(I)"
		     	#value = 0.001
		end
	push!(summarize,[mSet,"Random (Marker)",str,M[mSet][:df],M[mSet][:scale]])
	end
	

	push!(summarize,["e","Random",priorVCV[:e].str,dfE,scaleE])						

	println("\n ---------------- Summary of analysis ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,nosubheader=true,alignment=:l)


	#########make MCMC output files.
	
	levelsX = hcat([value[:levels] for (key, value) in X]...)
#	[println(key, value[:levels]) for (key, value) in X]
			
	IO.outMCMC(outPut,"b",levelsX)
	
	#check for correlated RE
        for zSet in keys(Z)
		if isa(zSet, Symbol)
			nameRE_VCV = String(zSet)
			IO.outMCMC(outPut,"u$zSet",[Z[zSet][:levels]])
			IO.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
		elseif isa(zSet, Expr)
			nameRE_VCV = join(zSet.args)[2:end]
			IO.outMCMC(outPut,"u$zSet",[Z[zSet][:levels]])
			IO.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
		elseif isa(zSet, Tuple)
			nameRE_VCV =  join(String.(vcat(zSet...)),"_").*hcat(["_$i" for i in 1:(length(zSet)^2)]...)
			for z in zSet
   				IO.outMCMC(outPut,"u$z",[Z[zSet][:levels]])
			end
			IO.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
		end
	end	
		
	#arbitrary marker names
	for mSet in keys(M)
		if isa(mSet,Tuple)
			for m in mSet
   				IO.outMCMC(outPut,"beta$m",hcat(M[mSet][:levels]...))
			end
		end
        end
	
	for mSet in keys(varBeta)
		isa(mSet, Symbol) ? nameM_VCV = ["reg_$r" for r in 1:M[mSet][:nRegions]] : nameM_VCV = vcat([["reg_$(i)_$j" for j in 1:size(M[mSet][:scale],2)^2] for i in 1:M[mSet][:nRegions]]...)
		IO.outMCMC(outPut,"var$mSet",[nameM_VCV]) #[] to have it as one row
        end
	

	IO.outMCMC(outPut,"varE",["e"])
	##########
	
	X  = myUnzip(X)
	Z  = myUnzip(Z)
	M  = myUnzip(M)	
	
	return ycorr, nData, dfE, scaleE, X, b, Z, u, varU, M,  beta, varBeta, BayesX
	
end

end
