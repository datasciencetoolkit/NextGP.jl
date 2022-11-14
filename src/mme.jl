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


#main sampler
function getMME!(iA,iGRel,Y,X,Z,M,levelDict,blocks,priorVCV,summaryStat,outPut)
		
        #some info
	nRand = length(Z)
	nColEachZ    = OrderedDict(z => size(Z[z],2) for z in keys(Z))
	nData = length(Y)
	nMarkerSets = length(M)
        #initial computations and settings
	ycorr = deepcopy(Y)

	priorVCV = convert(Dict{ExprOrSymbol, Any},priorVCV)	
	
	### X and b
	
	
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
	
	for blk in blocks
#		getThese = intersect(collect(keys(X)), blk)
		X[blk] = Dict{Symbol, Any}()
		X[blk][:data] = hcat(getindex.(getindex.(Ref(X), blk),:data)...)
		X[blk][:levels] = hcat(vcat(getindex.(getindex.(Ref(X), blk),:levels)...)...)
		X[blk][:nCol] = sum(getindex.(getindex.(Ref(X), blk),:nCol))
		X[blk][:pos] = vcat(getindex.(getindex.(Ref(X), blk),:pos)...)
		X[blk][:method] = first(getindex.(getindex.(Ref(X), blk),:method))
		for d in blk
			delete!(X,d)
		end
	end


#	b = zeros(getindex.(getindex.(Ref(X), keys(X)),:nCol)[])
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
					
        uKeyPos = OrderedDict{Any,Int64}()
        for zSet in keys(Z)
		pos = findall(x->x==zSet, collect(keys(Z)))[]
                uKeyPos[zSet] = pos
        end

	#matrices are ready
				
	Zp = OrderedDict{Any,Any}()
       	zpz = OrderedDict{Any,Any}() #Has the order in priorVCV, which may be unordered Dict() by the user. Analysis follow this order.
	rhsZ = OrderedDict{Any,Any}()
					
	for pSet ∈ keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV) 
		corEffects = []
		corPositions = []
		#symbol :ID or expression :(1|ID)
		if (isa(pSet,Symbol) || isa(pSet,Expr)) && in(pSet,keys(Z))
			tempzpz = []
			nowZ = Z[pSet]
			for c in eachcol(nowZ)
				push!(tempzpz,c'c)					
				# push!(tempzpz,BLAS.dot(c,c))
			end
			Zp[pSet]  = transpose(Z[pSet])						
			zpz[pSet] = tempzpz
			rhsZ[pSet] = zeros(size(Z[pSet],2))
                        if pSet in keys(summaryStat)
                                summaryStat[pSet].v == Array{Float64,1} ? zpz[pSet] .+= inv.(summaryStat[pSet].v) : zpz[pSet] .+= inv.(diag(summaryStat[pSet].v))
                                summaryStat[pSet].v == Array{Float64,1} ? rhsZ[pSet] .= inv.(summaryStat[pSet].v) .* (summaryStat[pSet].m)  : rhsZ[pSet] .= inv.(diag(summaryStat[pSet].v)) .* (summaryStat[pSet].m)
                        end
		#tuple of symbols (:ID,:Dam)
		elseif (isa(pSet,Tuple{Vararg{Symbol}})) && all((in).(pSet,Ref(keys(Z)))) #if all elements are available # all([pSet .in Ref(keys(Z))])
			correlate = collect(pSet)
			for pSubSet in correlate
				push!(corEffects,pSubSet)
				push!(corPositions,findall(pSubSet.==keys(Z))[])
			end
			if issubset(corEffects,collect(keys(Z)))
				tempZ = hcat.(eachcol.(getindex.(Ref(Z), (pSet)))...)
				for d in corEffects
                       			delete!(Z,d)
					delete!(uKeyPos,d)												
               			end
				uKeyPos[pSet] = corPositions
				Z[pSet]   = tempZ
				zpz[pSet] = MatByMat.(tempZ)
				Zp[pSet]  = transpose.(tempZ)
				tempZ = 0
				if pSet in keys(summaryStat)
					error("Not available to use summary statistics in correlated effects")
                        		#SummaryStat[pSet].v == Array{Float64,1} ? zpz[pSet] += inv.(SummaryStat[pSet].v) : zpz[pSet] += inv.(diag(SummaryStat[pSet].v))
                        		#SummaryStat[pSet].v == Array{Float64,1} ? rhsZ[pSet] = inv.(SummaryStat[pSet].v) .* (SummaryStat[pSet].m)  : rhsZ[pSet] = inv.(diag(SummaryStat[pSet].v)) .* (SummaryStat[pSet].m)
                		end
			end
		end
	end
																
	for pSet in collect(keys(Z))[(!in).(keys(Z),Ref(keys(priorVCV)))]
		printstyled("No prior was provided for $pSet, but it was not included in the data. It will be made uncorrelated with default priors\n"; color = :green)		
		tempzpz = []
		nowZ = Z[pSet]
		for c in eachcol(nowZ)
			push!(tempzpz,c'c)					
		end
		Zp[pSet]  = transpose(Z[pSet])						
		zpz[pSet] = tempzpz
		rhsZ[pSet] = zeros(size(Z[pSet],2))
		if pSet in keys(summaryStat)
                	summaryStat[pSet].v == Array{Float64,1} ? zpz[pSet] .+= inv.(summaryStat[pSet].v) : zpz[pSet] .+= inv.(diag(summaryStat[pSet].v))
                        summaryStat[pSet].v == Array{Float64,1} ? rhsZ[pSet] .= inv.(summaryStat[pSet].v) .* (summaryStat[pSet].m)  : rhsZ[pSet] .= inv.(diag(summaryStat[pSet].v)) .* (summaryStat[pSet].m)
                end
	end
																	
	#pos for individual random effect
	#this part "collect(k) .=> collect(v)" will change for correlated random effects.
	uKeyPos4Print = OrderedDict(vcat([(isa(k,Symbol) || isa(k,Expr)) ? k => v : collect(k) .=> collect(v) for (k,v) in uKeyPos]...))
	
	##get priors per effect
													
	iVarStr = Dict{Any,Array{Float64,2}}() #inverses will be computed
	varU_prior = OrderedDict{Any,Any}()
        for zSet in keys(Z)
                nCol = size(Z[zSet],2)
		#var structures and priors
		if haskey(priorVCV,zSet)	
			if ismissing(priorVCV[zSet].str) || priorVCV[zSet].str=="I" 
				printstyled("prior var-cov structure for $zSet is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				iVarStr[zSet] = Matrix(1.0I,nCol,nCol)
			elseif priorVCV[zSet].str=="A"
				printstyled("prior var-cov structure for $zSet is A. Computed A matrix (from pedigree file) will be used\n"; color = :green)
				iVarStr[zSet] = iA
			elseif priorVCV[zSet].str=="G"
                                printstyled("prior var-cov structure for $zSet is G. Computed G matrix will be used\n"; color = :green)
                                iVarStr[zSet] = iGRel[zSet]
			else 	iVarStr[zSet] = inv(priorVCV[zSet].str)
			end
			varU_prior[zSet] = priorVCV[zSet].v
		else	
			printstyled("prior var-cov for $zSet is empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
			varU_prior[zSet] = 100
			priorVCV[zSet] = Random("I",0,100)
			iVarStr[zSet] = Matrix(1.0I,nCol,nCol)
		end
        end

	#df, shape, scale...															
	
	dfZ = Dict{Any,Any}()	
	for zSet ∈ keys(zpz)
		dfZ[zSet] = 3.0+size(priorVCV[zSet].v,1)
	end
																
	scaleZ = Dict{Any,Any}()
        for zSet in keys(zpz)
                nZComp = size(priorVCV[zSet].v,1)
		#priorVCV[zSet].v is a temporary solution
		nZComp > 1 ? scaleZ[zSet] = priorVCV[zSet].v .* (dfZ[zSet]-nZComp-1.0)  : scaleZ[zSet] = priorVCV[zSet].v * (dfZ[zSet]-2.0)/dfZ[zSet] #I make float and array of float														
        end


												
        ####
																					

	#ADD MARKERS
	# read map file and make regions
																		
	############priorVCV cannot be empty for markers, currently!!																	

	#key positions for each effect in beta, for speed. Order of matrices in M are preserved here.
        for mSet in keys(M)
                pos = findall(mSet.==collect(keys(M)))[]
                M[mSet][:pos] = pos
        end

	beta = []

	#make mpm
	
	for pSet ∈ keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV)
		#symbol :M1 or expression
		if isa(pSet,Symbol) && in(pSet,keys(M))
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
			maps = getindex.(getindex.(Ref(M),pSet),:map)
			(length(maps)==0 || all( ==(maps[1]), maps)) == true ? M[pSet][:map] = first(maps) : error("correlated marker sets must have the same map file!")
			M[pSet][:pos] = vcat(getindex.(getindex.(Ref(M), pSet),:pos)...)
			M[pSet][:levels] = first(getindex.(getindex.(Ref(M),pSet),:levels))
			tempM = hcat.(eachcol.(getindex.(getindex.(Ref(M), pSet),:data))...)
			M[pSet][:data] = tempM
			for d in pSet
                       		delete!(M,d)
               		end
			M[pSet][:mpm] = MatByMat.(tempM)
			if pSet in keys(summaryStat)
				error("Not available to use summary statistics in correlated effects")
                                #SummaryStat[pSet].v == Array{Float64,1} ? mpm[pSet] += (1.0 ./ SummaryStat[pSet].v) : mpm[pSet] += inv.(diag(SummaryStat[pSet].v))
 	                end
			M[pSet][:Mp]  = transpose.(tempM)
			tempM = 0
			theseRegions = prep2RegionData(outPut,pSet,M[pSet][:map],priorVCV[pSet].r)
			M[pSet][:regionArray] = theseRegions
			M[pSet][:nRegions] = length(theseRegions)
		end
	end
	
	for pSet in collect(keys(M))[(!in).(keys(M),Ref(keys(priorVCV)))]
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
	u = zeros(Float64,nRand,maximum(vcat([0,collect(values(nColEachZ))]...))) #zero is for max to work when no random effect is present #can allow unequal length! Remove tail zeros for printing....

	varU = deepcopy(varU_prior) #for storage

	varBeta = Dict{Union{Symbol,Tuple{Vararg{Symbol}}},Any}()
        for mSet in keys(M)
                varBeta[mSet] = [priorVCV[mSet].v for i in 1:M[mSet][:nRegions]] #later, direct reference to key when varM_prior is a dictionary
        end

	#summarize analysis
	summarize = DataFrame(Effect=Any[],Type=Any[],Str=Any[],df=Any[],scale=Any[])
	
	for zSet in keys(zpz)
		if zSet ∈ keys(priorVCV)
			str = priorVCV[zSet].str
			#value = priorVCV[zSet].v
		else 
			str = "I"
		     	#value = varU_prior[zSet].v
		end
	push!(summarize,[zSet,"Random",str,dfZ[zSet],scaleZ[zSet]])
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
	#not a dictionary anymore, and consistent with possible new order.
	for (key, value) in X
		println("valueLevels: $(value[:levels])")
	end
	
	levelsX = hcat(vcat([value[:levels] for (key, value) in X]...)...)
	println("levelsX: $levelsX")
			
	IO.outMCMC(outPut,"b",levelsX)
	
	#check for correlated RE
        for i in 1:length(levelDict[:levelsRE])
		levRE = hcat(vcat(collect(values(levelDict[:levelsRE]))[i]...)...)
		IO.outMCMC(outPut,"u$i",levRE)
		isa(collect(keys(levelDict[:levelsRE]))[i], Symbol) ? nameRE_VCV = String(collect(keys(levelDict[:levelsRE]))[i]) : nameRE_VCV = join(collect(keys(levelDict[:levelsRE]))[i].args)[2:end]
		IO.outMCMC(outPut,"varU$i",[nameRE_VCV]) #[] to have it as one row
	end	
		
	#arbitrary marker names
	for mSet in keys(M)
   		IO.outMCMC(outPut,"beta$mSet",hcat(M[mSet][:levels]...))
        end
	
	for mSet in keys(varBeta)
		isa(mSet, Symbol) ? nameM_VCV = ["reg_$r" for r in 1:M[mSet][:nRegions]] : nameM_VCV = vcat([["reg_$(i)_$j" for j in 1:size(M[mSet][:scale],2)^2] for i in 1:M[mSet][:nRegions]]...)
		IO.outMCMC(outPut,"var$mSet",[nameM_VCV]) #[] to have it as one row
        end
	

	IO.outMCMC(outPut,"varE",["e"])
	##########
	
	X  = myUnzip(X)
	M  = myUnzip(M)	
	
	return ycorr, nData, dfE, scaleE, X, b, Z, iVarStr, Zp, zpz, uKeyPos, uKeyPos4Print, nColEachZ, u, varU, scaleZ, dfZ, M,  beta, varBeta, BayesX, rhsZ
	
end

end
