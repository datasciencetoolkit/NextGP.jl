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
ExprOrSymbol = Union{Expr,Symbol}
ExprOrSymbolOrTuple = Union{Expr,Symbol,Tuple}

#set up (co)variance structures
function setVarCovStr!(zSet::ExprOrSymbolOrTuple,Z::Dict,priorVCV,varU_prior::Dict)
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
		priorVCV[zSet] = Random("I",100)
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

	priorVCV = convert(Dict{ExprOrSymbolOrTuple, Any},priorVCV)	


	#set up for E.
	E = Dict{Any,Any}()	
	#no inverse implemented yet!
	if haskey(priorVCV,:e)	
		if isempty(priorVCV[:e].str) || priorVCV[:e].str=="I" 
				printstyled("prior var-cov structure for \"e\" is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				E[:str] = "I"
				E[:iVarStr] = [] #Matrix(1.0I,nData,nData)
				priorVCV[:e] = Random("I",priorVCV[:e].v)
		elseif isa(priorVCV[:e].str,Vector) # D
				E[:str] = "D"
				E[:iVarStr] = inv.(priorVCV[:e].str) #inv(Diagonal(priorVCV[:e].str))
#				error("var-cov structure \"D\" has not been implemented yet")
				printstyled("prior var-cov structure for \"e\" is \"D\". User provided \"D\" matrix (d_ii = 1/w_ii) will be used\n"; color = :green)
		else 
				error("provide a valid prior var-cov structure (\"I\", \"D\" or leave it empty \"[]\") for \"e\" ")
		end
	else	
		printstyled("prior var-cov for \"e\" is fully  empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		E[:iVarStr] = [] #Matrix(1.0I,nData,nData)
		#just add to priors
		priorVCV[:e] = Random("I",100)
	end
								
	#parameters for priors
        E[:df] = 4.0
 	       
	if priorVCV[:e].v==0.0
		priorVCV[:e].v  = 0.0005
       		E[:scale]     = 0.0005
        else
       		E[:scale]    = priorVCV[:e].v*(E[:df]-2.0)/E[:df]    
   	end
	
	### X and b	
	
	for blk in blocks
		println("blocking variables in $blk")
		X[blk] = Dict{Symbol, Any}()
		X[blk][:data] = hcat(getindex.(getindex.(Ref(X), blk),:data)...)
		X[blk][:levels] = vcat(getindex.(getindex.(Ref(X), blk),:levels)...)
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

	###Allow no fixed effects
	isempty(keys(X)) ? b = [] : b = zeros(sum(getindex.(getindex.(Ref(X), keys(X)),:nCol)))


	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	
        for xSet in keys(X)
		if E[:str] == "D"
#			X[xSet][:xpx] = X[xSet][:data]'*E[:iVarStr]*X[xSet][:data]
			X[xSet][:xpx] = X[xSet][:data]'*(E[:iVarStr].*X[xSet][:data])
			X[xSet][:Xp] = transpose(X[xSet][:data].*E[:iVarStr])
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
			setVarCovStr!(zSet,Z,priorVCV,varU_prior)
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
						else throw(ArgumentError("Please enter a valid region size (1 or 9999)"))
						end
					else
						theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],priorVCV[pSet].r)
						M[pSet][:regionArray] = theseRegions
					end
					M[pSet][:nVarCov] = length(theseRegions)
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
					designMat = modelmatrix(priorVCV[pSet].f, priorVCV[pSet].covariates)
					M[pSet][:covariates] = designMat
					M[pSet][:covariatesT] = transpose(designMat)
					M[pSet][:c] = rand(size(designMat,2))
					M[pSet][:SNPVARRESID] = rand(size(designMat,1))
					#iCpC inverse taken later
					M[pSet][:iCpC] = M[pSet][:covariatesT]*M[pSet][:covariates]
					if isa(M[pSet][:iCpC],Matrix{Float64}) 
						M[pSet][:iCpC] += Matrix(I*minimum(abs.(diag(M[pSet][:iCpC])./10000)),size(M[pSet][:iCpC]))
						#Matrix(I*0.001,size(M[pSet][:iCpC]))
					end
 		              		M[pSet][:iCpC]  = inv(M[pSet][:iCpC])
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
		
	for mSet ∈ keys(M)
		M[mSet][:df] = haskey(priorVCV,mSet) ? 3.0+size(priorVCV[mSet].v,1) : 4.0
		#M[mSet][:df] = 3.0+size(priorVCV[mSet].v,1)
	end


        for mSet ∈ keys(M)
		if haskey(priorVCV,mSet)
                	nMComp = size(priorVCV[mSet].v,1)
                	M[mSet][:scale] = nMComp>1 ? priorVCV[mSet].v .* (M[mSet][:df]-nMComp-1.0)  : priorVCV[mSet].v * (M[mSet][:df]-2.0)/(M[mSet][:df]) #I make float and array of float
		else
			nMComp = 1
			M[mSet][:scale] = 0.05 * (M[mSet][:df]-2.0)/(M[mSet][:df]) #I make float and array of floa
		end
        end
	
	
	#storage

	varU = deepcopy(varU_prior) #for storage

	varBeta = Dict{Union{Symbol,Tuple{Vararg{Symbol}}},Any}()
        for mSet ∈ keys(M)
		if haskey(priorVCV,mSet)
                	varBeta[mSet] = [priorVCV[mSet].v for i in 1:M[mSet][:nVarCov]]
		else
			varBeta[mSet] = [0.05 for i in 1:M[mSet][:nVarCov]]
		end
        end

	#summarize analysis
	summarize = DataFrame(Effect=Any[],Type=Any[],Str=Any[],df=Any[],scale=Any[])
	
	for zSet in keys(Z)
		push!(summarize,[zSet,"Random",Z[zSet][:str],Z[zSet][:df],Z[zSet][:scale]])
	end

	for mSet in keys(M)
		M[mSet][:method] == "BayesPR" ? str = "$(M[mSet][:method]) $(M[mSet][:nVarCov]) block(s)" : str = "$(M[mSet][:method])"
		push!(summarize,[mSet,"Random (Marker)",str,M[mSet][:df],M[mSet][:scale]])
	end
	

	push!(summarize,["e","Random",E[:str],E[:df],E[:scale]])						

	println("\n ---------------- Summary of analysis ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,nosubheader=true,alignment=:l)


	#########make MCMC output files.
	
	isempty(blocks) ? levelsX = hcat(vcat([value[:levels] for (key, value) in X]...)...) : levelsX = hcat(vcat([vcat(value[:levels]) for (key, value) in X]...)...)

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
		if isa(mSet,Symbol)
			IO.outMCMC(outPut,"beta$mSet",hcat(M[mSet][:levels]...))
			IO.outMCMC(outPut,"delta$mSet",hcat(M[mSet][:levels]...))
			if in(M[mSet][:method],["BayesB","BayesC","BayesR"])
				IO.outMCMC(outPut,"pi$mSet",[["pi$v" for v in 1:length(M[mSet][:vClass])]]) #[] to have it as one row
			elseif in(M[mSet][:method],["BayesRCπ","BayesRCplus"])
				npis = length(M[mSet][:vClass])*M[mSet][:nVarCov]
				IO.outMCMC(outPut,"pi$mSet",[["pi$v" for v in 1:npis]]) #[] to have it as one row
				IO.outMCMC(outPut,"annot$mSet",hcat(M[mSet][:levels]...))
			elseif in(M[mSet][:method],["BayesLV"])
				IO.outMCMC(outPut,"c$mSet",[["c$v" for v in 1:(length(M[mSet][:c]))]]) #[] to have it as one row
				IO.outMCMC(outPut,"varZeta$mSet",["varZeta"])
			end
		elseif isa(mSet,Tuple)
			for m in mSet
   				IO.outMCMC(outPut,"beta$m",hcat(M[mSet][:levels]...))
				IO.outMCMC(outPut,"delta$m",hcat(M[mSet][:levels]...))
			end
		end
        end
	
	for mSet in keys(varBeta)
		isa(mSet, Symbol) ? nameM_VCV = ["reg_$r" for r in 1:M[mSet][:nVarCov]] : nameM_VCV = vcat([["reg_$(i)_$j" for j in 1:size(M[mSet][:scale],2)^2] for i in 1:M[mSet][:nRegions]]...)
		IO.outMCMC(outPut,"var$mSet",[nameM_VCV]) #[] to have it as one row
        end
	

	IO.outMCMC(outPut,"varE",["e"])
	##########
	
	X  = myUnzip(X)
	Z  = myUnzip(Z)
	M  = myUnzip(M)
	E  = (;E...)
	
	return ycorr, nData, E, X, b, Z, u, varU, M,  beta, varBeta, delta
	
end

end
