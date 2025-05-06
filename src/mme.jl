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
include("types.jl")
include("varComp.jl")

#function name attached to genomic component, such as M[pSet][:funct] = sampleBayesC!
include("functions.jl")
using .functions

export getMME!


function MMEX!(X,eSet,E,blocks,summaryStat) #LHS is a Tuple
	println("dealing trait $eSet")
	if haskey(blocks, eSet)
		for blk in blocks[eSet]
			println("blocking variables $blk for trait $eSet")
			X[blk] = Dict{Symbol, Any}()
			X[blk][:data] = hcat(getindex.(getindex.(Ref(X), blk),:data)...)
			X[blk][:levels] = vcat(getindex.(getindex.(Ref(X), blk),:levels)...)
			X[blk][:nCol] = sum(getindex.(getindex.(Ref(X), blk),:nCol))
			X[blk][:method] = first(getindex.(getindex.(Ref(X), blk),:method))
			for d in blk
				delete!(X,d)
			end
		end
	else println("NO blocking is performed for trait $eSet")
	end
	

	#==BLOCK FIXED EFFECTS.
	Order of blocks is as defined by the user
	Order of variables within blocks is always the same as in the model definition, not defined by the user in each block.
	==#
	
        for xSet in keys(X)
		println("eSet: $eSet xSet: $xSet")
		if E[eSet][:str] == "D"
#			X[xSet][:xpx] = X[xSet][:data]'*E[ySet][:iVarStr]*X[xSet][:data]
			X[xSet][:xpx] = X[xSet][:data]'*(E[ySet][:iVarStr].*X[xSet][:data])
			X[xSet][:Xp] = transpose(X[xSet][:data].*E[ySet][:iVarStr])
		else 
			X[xSet][:xpx] = X[xSet][:data]'X[xSet][:data]
			X[xSet][:Xp] = transpose(X[xSet][:data])
		end
		X[xSet][:lhs] = zeros(X[xSet][:nCol])
		X[xSet][:rhs] = zeros(X[xSet][:nCol])

                if xSet in keys(summaryStat)
	  		X[xSet][:lhs] .= isa(summaryStat[xSet].v,Array{Float64,1}) ? inv.(summaryStat[xSet].v) : inv.(diag(summaryStat[xSet].v))
			X[xSet][:rhs] .= isa(summaryStat[xSet].v,Array{Float64,1}) ? inv.(summaryStat[xSet].v) .* (summaryStat[xSet].m)  : inv.(diag(summaryStat[xSet].v)) .* (summaryStat[xSet].m)
                end

		if isa(X[xSet][:xpx],Matrix{Float64})
#			println("diag: $(diag(X[xSet][:xpx])) added to diag: $(minimum(abs.(diag(X[xSet][:xpx]))))")
			X[xSet][:xpx] += Matrix(I*minimum(abs.(diag(X[xSet][:xpx])./10000)),size(X[xSet][:xpx]))
		end
	end
end



#main sampler
function getMME!(Y,X,Z,M,E,blocks,priorVCV,summaryStat,modelInformation,outPut) #maybe later use modelInformation
		
        #some info
	nRand = length(Z)
	nData = length(Y)
	nMarkerSets = length(M)
        
	#initial computations and settings
	
	ycorr = deepcopy(Y)

	priorVCV = convert(Dict{ExprOrSymbolOrTuple, Any},priorVCV)
	
	
	varU = Dict{Any,Any}() #for storage
	varBeta = Dict{Union{Symbol,Tuple{Vararg{Symbol}}},Any}()
	varE = Dict{Union{Symbol,Tuple{Vararg{Symbol}}},Any}()

	######## 
	#E (is per trait information up until now)
	########


	#set up varCov for e
	for eSet in keys(E)
		setVarCovStrE!(eSet,E,priorVCV,nData,varE)
	end
	varCovE!(E,priorVCV)

	println(E)

	###################################
	
	
	######## 
	#X and b
	########
	
	if isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Symbol
		println("model is a single-trait model")
		for eSet in keys(E)
			MMEX!(X,eSet,E,blocks,summaryStat)
		end
	elseif isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Tuple
		println("model is a multi-trait model where measurements/observations are from the same individuals")
	elseif !isequal(length(collect(keys(E))),1) && all(typeof.(collect(keys(E))) .<: Symbol)
		println("model is a multi-population model where measurements/observations are from different individuals")
	else 	ArgumentError("Could not understand the type of your model")

	end

	#Positions of parameters for each variable and blocks for speed. b is a column vector.
	countXCol = 0
	for xSet in keys(X)
		nCol = X[xSet][:nCol]
          	X[xSet][:pos] = (countXCol+1):(countXCol+nCol)
		countXCol += nCol
        end
	countXCol = 0

	###Allow no fixed effects
	isempty(keys(X)) ? b = [] : b = zeros(sum(getindex.(getindex.(Ref(X), keys(X)),:nCol)))
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)

	###################################
	
	### 
	#Z and u
	###
	
	u = []

	#matrices are ready
	
	posZcounter = 0
	for zSet in keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV)
		#symbol :ID or expression :(1|ID)
		if (isa(zSet,Symbol) || isa(zSet,Expr)) && in(zSet,keys(Z))
			posZcounter += 1
			Z[zSet][:pos] = posZcounter
			tempzpz = []
			nowZ = Z[zSet][:data]
			
			setVarCovStrU!(zSet,Z,priorVCV,varU_prior)
			
			if E[:str] == "D"
				for c in eachcol(nowZ)
#					push!(tempzpz,dot(c,E[:iVarStr],c))
					push!(tempzpz,sum(c.*E[:iVarStr].*c))
				end
#				Z[zSet][:Zp]  = transpose(nowZ)*E[:iVarStr]
				Z[zSet][:Zp]  = transpose(nowZ.*E[:iVarStr])
			else
				for c in eachcol(nowZ)
					push!(tempzpz,dot(c,c))
				end
				Z[zSet][:Zp]  = transpose(nowZ)
			end
			Z[zSet][:zpz] = tempzpz
			u = push!(u,zeros(Float64,1,size(nowZ,2)))
			nowZ = 0
			tempzpz = 0

			#summary statistics
			Z[zSet][:lhs] = zeros(size(nowZ,2))
			Z[zSet][:rhs] = zeros(size(nowZ,2))
                        if zSet in keys(summaryStat)
                                Z[zSet][:lhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) : inv.(diag(summaryStat[zSet].v))
                                Z[zSet][:rhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) .* (summaryStat[zSet].m)  : inv.(diag(summaryStat[zSet].v)) .* (summaryStat[zSet].m)
                        end
			
		#tuple of symbols (:ID,:Dam)
		elseif (isa(zSet,Tuple{Vararg{Symbol}})) && all((in).(zSet,Ref(keys(Z)))) #if all elements are available # all([zSet .in Ref(keys(Z))])
			Z[zSet] = Dict{Symbol, Any}()
#			Z[zSet][:pos] = collect((posZcounter+1):(posZcounter+length(zSet)))
#			posZcounter += length(zSet)
			posZcounter += 1
			Z[zSet][:pos] = posZcounter
			u = push!(u,zeros(Float64,length(zSet),size(Z[zSet[1]][:data],2)))
			Z[zSet][:levels] = first(getindex.(getindex.(Ref(Z),zSet),:levels))
			tempZ = hcat.(eachcol.(getindex.(getindex.(Ref(Z), zSet),:data))...)
			#same Z for all components in a single-trait model get only first column! Z[zSet][:data] = getindex.(tempZ,:,1)
			Z[zSet][:data] = tempZ
			
			setVarCovStrU!(zSet,Z,priorVCV,varU_prior)
			
			Z[zSet][:str] = Z[zSet[1]][:str] 
			
			for d in zSet
                       		delete!(Z,d)
               		end

			###WEIGHTED SHOULD BE ADAAPTED HERE#################
			Z[zSet][:zpz]  = MatByMat.(tempZ)
			Z[zSet][:Zp]   = transpose.(tempZ)


			
			#lhs is already zero as only mpm + "nothing" is  given
			#rhs is for now only for convenience
			Z[zSet][:rhs] = [zeros(length(zSet)) for i in 1:length(Z[zSet][:levels])]
			if zSet in keys(summaryStat)
				error("Not available to use summary statistics in correlated effects")
                        	#SummaryStat[pSet].v == Array{Float64,1} ? zpz[pSet] += inv.(SummaryStat[pSet].v) : zpz[pSet] += inv.(diag(SummaryStat[pSet].v))
                        	#SummaryStat[pSet].v == Array{Float64,1} ? rhsZ[pSet] = inv.(SummaryStat[pSet].v) .* (SummaryStat[pSet].m)  : rhsZ[pSet] = inv.(diag(SummaryStat[pSet].v)) .* (SummaryStat[pSet].m)
                	end
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
		Z[zSet][:Zp]  = transpose(nowZ)						
		Z[zSet][:zpz] = tempzpz
		Z[zSet][:lhs] = zeros(size(nowZ,2))
		Z[zSet][:rhs] = zeros(size(nowZ,2))
		if zSet in keys(summaryStat)
                	Z[zSet][:lhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) : inv.(diag(summaryStat[zSet].v))
                        Z[zSet][:rhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) .* (summaryStat[zSet].m)  : inv.(diag(summaryStat[zSet].v)) .* (summaryStat[zSet].m)
                end
		setVarCovStr!(zSet,Z,priorVCV,varU_prior)
	end

	##set up varCov for u
	##varCovZ!(Z,priorVCV,varU_prior,varU)
																		
												
        ####
																					

	#ADD MARKERS
	# read map file and make regions
																		
	beta = []
	delta = []

	#make mpm
	posMcounter = 0
	for pSet in keys(M) 		#keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV)
#		in(pSet, collect(keys(M))[(!in).(keys(M),Ref(keys(priorVCV)))]) ? throw(ArgumentError("You must provide a prior for genomic analysis. Example: BayesPR(9999,0.05)")) : nothing
		
		haskey(priorVCV,pSet) ? nothing : printstyled("No prior was provided for $pSet, but it was included in the data. It will be made uncorrelated with default priors and region size 9999 (WG)\n"; color = :green)
		
		#symbol :M1 or expression
		if isa(pSet,Symbol) #&& in(pSet,keys(M)), now keys are always in keys(M)
			posMcounter += 1
			M[pSet][:pos] = posMcounter
			tempmpm = []
			nowM = M[pSet][:data]

			if E[:str] == "D"
				for c in eachcol(nowM)
					push!(tempmpm,sum(c.*E[:iVarStr].*c))
				end
				M[pSet][:Mp] = map(i -> transpose(nowM[:,i].*E[:iVarStr]), axes(nowM, 2))
			else
				for c in eachcol(nowM)
					push!(tempmpm,dot(c,c))
				end
				M[pSet][:Mp] = map(i -> transpose(nowM[:,i]), axes(nowM, 2))
			end			

			M[pSet][:mpm] = tempmpm
			
			#summary statistics
			M[pSet][:lhs] = zeros(M[pSet][:dims][2])
			M[pSet][:rhs] = zeros(M[pSet][:dims][2])
			if pSet in keys(summaryStat)
                       		M[pSet][:lhs] .= isa(summaryStat[pSet].v,Array{Float64,1}) ? inv.(summaryStat[pSet].v) : inv.(diag(summaryStat[pSet].v))
				M[pSet][:rhs] .= isa(summaryStat[pSet].v,Array{Float64,1}) ? inv.(summaryStat[pSet].v) .* (summaryStat[pSet].m)  : inv.(diag(summaryStat[pSet].v)) .* (summaryStat[pSet].m)
				####Deal with N(0,0)
				M[pSet][:lhs][isinf.(M[pSet][:lhs])].= 0.0
				M[pSet][:rhs][isnan.(M[pSet][:rhs])].= 0.0
			end

			if !haskey(priorVCV,pSet)
				M[pSet][:method] = "BayesPR"
				M[pSet][:funct] = sampleBayesPR!
				theseRegions = [1:r for r in size(nowM,2)]
				M[pSet][:regionArray] = theseRegions
				M[pSet][:nVarCov] = length(theseRegions)
			else
				if priorVCV[pSet].name == "BayesPR"
					M[pSet][:method] = "BayesPR"
					M[pSet][:funct] = sampleBayesPR!
					if isempty(M[pSet][:map])		
						if priorVCV[pSet].r == 1
							printstyled("No map was provided. Running Bayesian Random Regression (BRR) with 1 SNP region size\n"; color = :green)
							theseRegions = [r:r for r in 1:size(nowM,2)]
							M[pSet][:regionArray] = theseRegions
						elseif priorVCV[pSet].r == 9999
							printstyled("No map was provided. Running Bayesian Random Regression (BRR) with all SNP as 1 region\n"; color = :green)
							theseRegions = [1:r for r in size(nowM,2)]
							M[pSet][:regionArray] = theseRegions
						else throw(ArgumentError("Please enter a valid region size (1 or 9999) or provide a map file"))
						end
					else
						theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],priorVCV[pSet].r)
						M[pSet][:regionArray] = theseRegions
					end
					M[pSet][:nVarCov] = length(theseRegions)
					M[pSet][:scale]   = [] 
				elseif priorVCV[pSet].name == "BayesB"
					M[pSet][:logPi]       = [log(1.0 .- priorVCV[pSet].pi) log(priorVCV[pSet].pi)] #not fitted, fitted
#					M[pSet][:logPiIn]     = log(priorVCV[pSet].pi)
#					M[pSet][:logPiOut]    = log(1.0 .- priorVCV[pSet].pi)
					M[pSet][:method]      = "BayesB"
					M[pSet][:funct]       = sampleBayesB!
					theseRegions          = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
					M[pSet][:nVarCov]     = length(theseRegions)
					M[pSet][:estPi]       = priorVCV[pSet].estimatePi
					M[pSet][:piHat]       = [1.0 .- priorVCV[pSet].pi priorVCV[pSet].pi] #not fitted, fitted
					M[pSet][:vClass]      = [0 1] #2 variance class, one with own, one with null
				elseif priorVCV[pSet].name == "BayesC"
					M[pSet][:logPi]       = [log(1.0 .- priorVCV[pSet].pi) log(priorVCV[pSet].pi)] #not fitted, fitted
#					M[pSet][:logPiIn]     = log(priorVCV[pSet].pi)
#					M[pSet][:logPiOut]    = log(1.0 .- priorVCV[pSet].pi)
					M[pSet][:method]      = "BayesC"
					M[pSet][:funct] = sampleBayesC!
					theseRegions          = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
					M[pSet][:nVarCov]     = 1
					M[pSet][:estPi]       = priorVCV[pSet].estimatePi
					M[pSet][:piHat]       = [1.0 .- priorVCV[pSet].pi priorVCV[pSet].pi] #not fitted, fitted
					M[pSet][:vClass]      = [0 1] #2 variance class, one with common, one with null
				elseif priorVCV[pSet].name == "BayesR"
					M[pSet][:logPi]       = log.(priorVCV[pSet].pi)
					M[pSet][:vClass]      = priorVCV[pSet].class
					M[pSet][:method]      = "BayesR"
					M[pSet][:funct]       = sampleBayesR!
					theseRegions          = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
					M[pSet][:nVarCov]     = 1
					M[pSet][:estPi]       = priorVCV[pSet].estimatePi
					M[pSet][:piHat]       = deepcopy(priorVCV[pSet].pi)
				elseif priorVCV[pSet].name == "BayesRCπ"
					M[pSet][:vClass]      = priorVCV[pSet].class
					M[pSet][:method]      = "BayesRCπ"
					M[pSet][:funct]       = sampleBayesRCπ!
					theseRegions          = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
					M[pSet][:nVarCov]     = size(priorVCV[pSet].annot,2)
					M[pSet][:logPi]       = [log.(priorVCV[pSet].pi) for i in 1:M[pSet][:nVarCov]]
					M[pSet][:estPi]       = priorVCV[pSet].estimatePi
					M[pSet][:piHat]       = [priorVCV[pSet].pi for i in 1:M[pSet][:nVarCov]]
					M[pSet][:annotInput]  = deepcopy(priorVCV[pSet].annot)
					M[pSet][:annotProb]   = priorVCV[pSet].annot./sum(priorVCV[pSet].annot,dims=2)
					#If all annotations are zero, prob is NA. I make it "0" here
					#But if all zero, those SNPs should be in a seperate marker set.
					#So i cancel this
					#M[pSet][:annotProb][findall([all(iszero, row) for row in eachrow(priorVCV[pSet].annot)]),:] .= 0.0
					#
					M[pSet][:annotNonZeroPos]   = [findall(!iszero, row) for row in eachrow(priorVCV[pSet].annot)]
#					M[pSet][:annotNonZero]= getindex.(Ref(priorVCV[pSet].annot),M[pSet][:annotNonZeroPos])
					M[pSet][:annotCat]    = zeros(Int64,1,M[pSet][:dims][2])
				elseif priorVCV[pSet].name == "BayesRCplus"
					M[pSet][:vClass]      = priorVCV[pSet].class
					M[pSet][:method]      = "BayesRCplus"
					M[pSet][:funct]       = sampleBayesRCplus!
					theseRegions          = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
					M[pSet][:nVarCov]     = size(priorVCV[pSet].annot,2)
					M[pSet][:logPi]       = [log.(priorVCV[pSet].pi) for i in 1:M[pSet][:nVarCov]]
					M[pSet][:estPi]       = priorVCV[pSet].estimatePi
					M[pSet][:piHat]       = [priorVCV[pSet].pi for i in 1:M[pSet][:nVarCov]]
					M[pSet][:annotInput]  = deepcopy(priorVCV[pSet].annot)
					M[pSet][:annotProb]   = priorVCV[pSet].annot./sum(priorVCV[pSet].annot,dims=2)
					M[pSet][:annotNonZeroPos]   = [findall(!iszero, row) for row in eachrow(priorVCV[pSet].annot)]
					M[pSet][:annotCat]    = zeros(Int64,1,M[pSet][:dims][2])
				elseif priorVCV[pSet].name == "BayesLV"
					M[pSet][:method]      = "BayesLV"
					M[pSet][:funct]       = sampleBayesLV!
					theseRegions          = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
					M[pSet][:nVarCov]     = length(theseRegions)
					#logVar can be created in a smarter way, maybe together with var??...
					M[pSet][:logVar]     = [log(priorVCV[pSet].v) for i in 1:M[pSet][:nVarCov]]
				#	designMat = modelmatrix(priorVCV[pSet].f, priorVCV[pSet].covariates)
				#	M[pSet][:covariates] = designMat
				#	M[pSet][:covariatesT] = transpose(designMat)
				#	M[pSet][:c] = rand(size(designMat,2))
				#	M[pSet][:SNPVARRESID] = rand(size(designMat,1))
					#iCpC inverse taken later
				#	M[pSet][:iCpC] = M[pSet][:covariatesT]*M[pSet][:covariates]
				#	if isa(M[pSet][:iCpC],Matrix{Float64}) 
				#		M[pSet][:iCpC] += Matrix(I*minimum(abs.(diag(M[pSet][:iCpC])./10000)),size(M[pSet][:iCpC]))
						#Matrix(I*0.001,size(M[pSet][:iCpC]))
				#	end
 		              	#	M[pSet][:iCpC]  = inv(M[pSet][:iCpC])
					M[pSet][:varZeta]  = [priorVCV[pSet].varZeta]
					M[pSet][:estVarZeta] = priorVCV[pSet].estimateVarZeta
					designMat = 0
				end
			end
			beta  = push!(beta,zeros(Float64,1,M[pSet][:dims][2]))
			delta = push!(delta,ones(Int64,1,M[pSet][:dims][2]))
			tempmpm = 0
			nowM = 0
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
																
			if isempty(M[pSet][:map])		
				if priorVCV[pSet].r == 1
					printstyled("No map was provided. Running Bayesian Random Regression (BRR) with 1 SNP region size\n"; color = :green)
					theseRegions = [r:r for r in 1:size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
				elseif priorVCV[pSet].r == 9999
					printstyled("No map was provided. Running Bayesian Random Regression (BRR) with all SNP as 1 region\n"; color = :green)
					theseRegions = [1:r for r in size(nowM,2)]
					M[pSet][:regionArray] = theseRegions
				else throw(ArgumentError("Please enter a valid region size (1 or 9999)"))
				end
			else
				theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],priorVCV[pSet].r)
				M[pSet][:regionArray] = theseRegions
			end
			M[pSet][:nVarCov] = length(theseRegions)
		end
	end

	##set up varCov for markers
	varCovM!(M,priorVCV,varBeta)	

	#summarize analysis
	summarize = DataFrame(Effect=Any[],Type=Any[],Str=Any[],df=Any[],scale=Any[])
	
	for zSet in keys(Z)
		push!(summarize,[zSet,"Random",Z[zSet][:str],Z[zSet][:df],Z[zSet][:scale]])
	end

	for mSet in keys(M)
		M[mSet][:method] == "BayesPR" ? str = "$(M[mSet][:method]) $(M[mSet][:nVarCov]) block(s)" : str = "$(M[mSet][:method])"
		push!(summarize,[mSet,"Random (Marker)",str,M[mSet][:df],M[mSet][:scale]])
	end
	
	for eSet in keys(E)
		push!(summarize,[eSet,"Random",E[eSet][:str],E[eSet][:df],E[eSet][:scale]])						
	end
	println("\n ---------------- Summary of analysis ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,alignment=:l)


	#########make MCMC output files.
	
	isempty(blocks) ? levelsX = hcat(vcat([value[:levels] for (key, value) in X]...)...) : levelsX = hcat(vcat([vcat(value[:levels]) for (key, value) in X]...)...)

	inOut.outMCMC(outPut,"b",levelsX)
	
	#check for correlated RE
        for zSet in keys(Z)
		if isa(zSet, Symbol)
			nameRE_VCV = String(zSet)
			inOut.outMCMC(outPut,"u$zSet",[Z[zSet][:levels]])
			inOut.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
		elseif isa(zSet, Expr)
			nameRE_VCV = join(zSet.args)[2:end]
			inOut.outMCMC(outPut,"u$zSet",[Z[zSet][:levels]])
			inOut.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
		elseif isa(zSet, Tuple)
			nameRE_VCV =  join(String.(vcat(zSet...)),"_").*hcat(["_$i" for i in 1:(length(zSet)^2)]...)
			for z in zSet
   				inOut.outMCMC(outPut,"u$z",[Z[zSet][:levels]])
			end
			inOut.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
		end
	end	
		
	#arbitrary marker names
	for mSet in keys(M)
		if isa(mSet,Symbol)
			inOut.outMCMC(outPut,"beta$mSet",hcat(M[mSet][:levels]...))
			inOut.outMCMC(outPut,"delta$mSet",hcat(M[mSet][:levels]...))
			if in(M[mSet][:method],["BayesB","BayesC","BayesR"])
				inOut.outMCMC(outPut,"pi$mSet",[["pi$v" for v in 1:length(M[mSet][:vClass])]]) #[] to have it as one row
			elseif in(M[mSet][:method],["BayesRCπ","BayesRCplus"])
				npis = length(M[mSet][:vClass])*M[mSet][:nVarCov]
				inOut.outMCMC(outPut,"pi$mSet",[["pi$v" for v in 1:npis]]) #[] to have it as one row
				inOut.outMCMC(outPut,"annot$mSet",hcat(M[mSet][:levels]...))
			elseif in(M[mSet][:method],["BayesLV"])
				inOut.outMCMC(outPut,"c$mSet",[["c$v" for v in 1:(length(M[mSet][:c]))]]) #[] to have it as one row
				inOut.outMCMC(outPut,"varZeta$mSet",["varZeta"])
			end
		elseif isa(mSet,Tuple)
			for m in mSet
   				inOut.outMCMC(outPut,"beta$m",hcat(M[mSet][:levels]...))
				inOut.outMCMC(outPut,"delta$m",hcat(M[mSet][:levels]...))
			end
		end
        end
	
	for mSet in keys(varBeta)
		isa(mSet, Symbol) ? nameM_VCV = ["reg_$r" for r in 1:M[mSet][:nVarCov]] : nameM_VCV = vcat([["reg_$(i)_$j" for j in 1:size(M[mSet][:scale],2)^2] for i in 1:M[mSet][:nRegions]]...)
		inOut.outMCMC(outPut,"var$mSet",[nameM_VCV]) #[] to have it as one row
		####
		isa(mSet, Symbol) ? nameM_VCV = "scale" : nameM_VCV = ["scale$i" for i in 1:size(vec(varBeta[mSet]))]
		inOut.outMCMC(outPut,"scale$mSet",[nameM_VCV]) #[] to have it as one row
		####
        end
	

	inOut.outMCMC(outPut,"varE",["e"])
	##########
	
	X  = myUnzip(X)
	Z  = myUnzip(Z)
	M  = myUnzip(M)
	E  = myUnzip(E) #(;E...)
	
	return ycorr, nData, E, varE, X, b, Z, u, varU, M,  beta, varBeta, delta
	
end

end
