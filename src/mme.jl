

#include("outFiles.jl")
#include("misc.jl")
#include("types.jl")

include("varComp.jl")
include("stBWGR.jl")

#function name attached to genomic component, such as M[mSet][:funct] = sampleBayesC!
#include("functions.jl")

#using .functions

#export getMME!


function blockX!(X,eSet,blocks,modelInformation) #LHS is a Tuple
	#println("modelInformation $modelInformation")
	#println("dealing trait $eSet")
	if haskey(blocks, eSet)
		#println("blocking variables $blocks for trait $eSet")
		for blk in blocks[eSet]
			#println("blocking variable $blk for trait $eSet")
			X[blk] = Dict{Symbol, Any}()
			X[blk][:data] = hcat(getindex.(getindex.(Ref(X), blk),:data)...)
			X[blk][:levels] = vcat(getindex.(getindex.(Ref(X), blk),:levels)...)
			X[blk][:nCol] = sum(getindex.(getindex.(Ref(X), blk),:nCol))
			X[blk][:method] = first(getindex.(getindex.(Ref(X), blk),:method))
			modelInformation[eSet][blk] = BlockTerm(blk)  #push!(collect(values(modelInformation[eSet])),blk)
			for d in blk
				println("deleting $d")
				delete!(X,d)
				delete!(modelInformation[eSet],d)
			end
			println("modelInformation $modelInformation FINAL")
		end
	else println("NO blocking is performed for trait $eSet")
	end
end
	

	#==BLOCK FIXED EFFECTS.
	Order of blocks is as defined by the user
	Order of variables within blocks is always the same as in the model definition, not defined by the user in each block.
	==#

#single-trait
function MMEX!(X,b,posXcounter,eSet::Symbol,E,blocks,modelInformation,summaryStat)
	#println("eSet is a Symbol")
	#println("X: $X")
	blockX!(X,eSet,blocks,modelInformation)
        for xSet in keys(X)
		posXcounter += 1 #position of this XSet's vector of effects in the big b vector
		X[xSet][:pos] = posXcounter
		if E[eSet][:str] == "D"
			X[xSet][:xpx] = X[xSet][:data]'*(E[ySet][:iVarStr].*X[xSet][:data])
			X[xSet][:Xp] = transpose(X[xSet][:data].*E[ySet][:iVarStr])
			#X[xSet][:Xp] = map(i -> transpose(X[xSet][:data][:,i].*E[eSet][:iVarStr]), axes(X[xSet][:data], 2))
		else 
			X[xSet][:xpx] = X[xSet][:data]'X[xSet][:data]
			X[xSet][:Xp] = transpose(X[xSet][:data])
			#X[xSet][:Xp] = map(i -> transpose(X[xSet][:data][:,i]), axes(X[xSet][:data], 2))
		end

		#summary statistics
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
		push!(b,zeros(Float64,X[xSet][:nCol],1))
	end
end


#multi-trait
function MMEX!(X,b,posXcounter,eSet::Tuple,E,blocks,modelInformation,summaryStat)
	#println("eSet: $eSet is a Tuple")
	blockX!(X,eSet,blocks,modelInformation)
	#println("X: $X")
        for xSet in keys(X)
		#println("MMEX: $xSet")
		posXcounter += 1 #position of this XSet's vector of effects in the big b vector
		X[xSet][:pos] = posXcounter
		
		xCol2Repeat = ntuple(i->xSet,length(eSet))
		#println("xCol2Repeat: $xCol2Repeat")

		tempX = hcat.(eachcol.(getindex.(getindex.(Ref(X), xCol2Repeat),:data))...)
		X[xSet][:data] = tempX #array of arrays where each array is a column in X[:xSet] for two traits size n*2 each of the arrays in the arrays]
		push!(b,[zeros(Float64,1,length(eSet)) for j in 1:X[xSet][:nCol]]) #is an array of arrays also!
						
		#Matrix of matrixces	
		X[xSet][:XpX] = hcat([[x'*tempX[j] for j in 1:length(tempX)] for x in tempX]...) #returns a big matrix of k'k!
		#println("size XpX: $xSet $(size(X[xSet][:XpX]))")
		
		#X[xSet][:XpX] = [[x'*tempX[j] for j in 1:length(tempX)] for x in tempX] #returns k array of arrays of t*t!
		#println("X[xSet][:XpX]: $(X[xSet][:XpX])")
		#X[xSet][:XpX] = [reduce(hcat, X[xSet][:XpX][i,:]) for i in 1:X[xSet][:nCol]]
		#X[xSet][:Xp]  = transpose.(tempX)
		#println("X[xSet][:XpX]: $(X[xSet][:XpX])")

		tempX = 0
		
		#if E[eSet][:str] == "D"
		#	for c in eachcol(nowX)
		#		#I compute x'X so that for each effect (column), I have both xpx and xpX stored! 
		#		push!(tempxpx,ones(length(eSet),length(eSet)).*sum(c.*E[eSet][:iVarStr].*c))
		#		push!(tempxpX,sum(c.*E[eSet][:iVarStr].*nowX,dims=1))
		#	end
		#	X[xSet][:Xp] = ones(length(eSet),1) .* map(i -> transpose(nowX[:,i].*E[eSet][:iVarStr]), axes(nowX, 2))
		#else
		#	for c in eachcol(nowX)
		#		#tempM = hcat.(eachcol.(getindex.(getindex.(Ref(X), xSet),:data))...)
		#		#I compute x'X so that for each effect (column), I have both xpx and xpX stored! 
		#		push!(tempxpx,ones(length(eSet),length(eSet)).*sum(c.*c))
		#		push!(tempxpX,sum(c.*nowX,dims=1))
		#	end
		#	X[xSet][:Xp] = [map(i -> transpose(nowX[:,i]), axes(nowX, 2)) for col in 1:length(eSet)]
		#end
		
		
		
		#summary statistics
		
		[X[xSet][:XpX][x] += Matrix(I*minimum(abs.(diag(X[xSet][:XpX][x])./10000)),size(X[xSet][:XpX][x])) for x in CartesianIndices(X[xSet][:XpX]) if isa(X[xSet][:XpX][x],Matrix{Float64})]
	end
end

#single-trait zSet::Union{Symbol,Tuple{Vararg{Symbol}}} allows for correlated effects 
#for multi-trait it will be, zSet::Union{Tuple{Tuple{Vararg{Symbol}}} #check potential overlap with single-trait
function MMEZ!(Z,u,varU,posZcounter,eSet::Symbol,E,priorVCV,modelInformation,summaryStat)
	#more like blockX function, but a bit different as the way it forms data structures.
	correlate = hcat(filter!(!isempty, unique([[keya for keya in keys(priorVCV) if (isa(keya,Tuple) && in(keyz,keya))] for keyz in keys(Z)]))...)
	#need to implement here data preparation w/o correlation
	if isempty(correlate)
		println("No Correlated Random effects")
	else 
		for zSet in correlate
			println("Correlating $correlate for $eSet")
			Z[zSet] = Dict{Symbol, Any}() #now Z has Dict(s) for the correlated effects
			Z[zSet][:iVarStr] = Z[zSet[1]][:iVarStr]
			modelInformation[eSet][zSet] = CorrelatedPedigreeTerm(zSet)
		end
	end
	
	for zSet in keys(Z)
		if (isa(zSet,Symbol) || isa(zSet,Expr)) #symbol :ID or expression :(1|ID)
			if any(in.(zSet,correlate))
				continue
			end
			posZcounter += 1 #should be here to aovid correlated ones having a positive first than be removed
			Z[zSet][:pos] = posZcounter
			
			tempzpz = []
			nowZ = Z[zSet][:data]
			if E[eSet][:str] == "D"
				for c in eachcol(nowZ)
					push!(tempzpz,sum(c.*E[eSet][:iVarStr].*c))
				end
				Z[zSet][:Zp]  = transpose(nowZ.*E[eSet][:iVarStr])
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
		elseif isa(zSet,Tuple)
			posZcounter += 1
			Z[zSet][:pos] = posZcounter
			Z[zSet][:levels] = first(getindex.(getindex.(Ref(Z),zSet),:levels))
			tempZ = hcat.(eachcol.(getindex.(getindex.(Ref(Z), zSet),:data))...)
			#same Z for all components in a single-trait model get only first column! Z[zSet][:data] = getindex.(tempZ,:,1)
			Z[zSet][:data] = tempZ
			Z[zSet][:str] = Z[zSet[1]][:str] 
			for d in zSet
                    	   	delete!(Z,d)
				delete!(modelInformation[eSet],d)
               		end
			u = push!(u,zeros(Float64,length(zSet),length(Z[zSet][:data])))
			###WEIGHTED SHOULD BE ADAAPTED HERE#################
			Z[zSet][:zpz]  = MatByMat.(tempZ)
			Z[zSet][:Zp]   = transpose.(tempZ)
			tempZ = 0
		else throw(ArgumentError("Could not understand the type of $zSet in Z"))
			
		end
	end

	for zSet in keys(Z)
		setVarCovStrU!(eSet,zSet,Z,priorVCV,varU)
	end
	
#	for zSet in collect(keys(Z))[(!in).(keys(Z),Ref(keys(priorVCV)))]
#		posZcounter += 1
#		Z[zSet][:pos] = posZcounter
#		printstyled("No prior was provided for $mSet, but it was not included in the data. It will be made uncorrelated with default priors\n"; color = :green)		
#		tempzpz = []
#		nowZ = Z[zSet][:data]
#		for c in eachcol(nowZ)
#			push!(tempzpz,c'c)					
#		end
#		Z[zSet][:Zp]  = transpose(nowZ)						
#		Z[zSet][:zpz] = tempzpz
#		Z[zSet][:lhs] = zeros(size(nowZ,2))
#		Z[zSet][:rhs] = zeros(size(nowZ,2))
#		if zSet in keys(summaryStat)
#                	Z[zSet][:lhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) : inv.(diag(summaryStat[zSet].v))
#                        Z[zSet][:rhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) .* (summaryStat[zSet].m)  : inv.(diag(summaryStat[zSet].v)) .* (summaryStat[zSet].m)
#                end
#	end
end

#single-trait, allows for correlated effects
function MMEM!(M,beta,varBeta,posMcounter,eSet::Symbol,E,priorVCV,modelInformation,summaryStat)
	#more like blockX function, but a bit different as the way it forms data structures.
	correlate = hcat(filter!(!isempty, unique([[keya for keya in keys(priorVCV) if (isa(keya,Tuple) && in(keyz,keya))] for keyz in keys(M)]))...)
	#need to implement here data preparation w/o correlation
	if isempty(correlate)
		println("No Correlated Marker effects")
	else 
		for mSet in correlate
			println("Correlating $correlate for $eSet")
			M[mSet] = Dict{Symbol, Any}() #now M has Dict(s) for the correlated effects
			#M[mSet][:iVarStr] = M[mSet[1]][:iVarStr]
			modelInformation[eSet][mSet] = CorrelatedMarkerTerm(mSet)
		end
	end
	
	for mSet in keys(M)
		if (isa(mSet,Symbol))
			if any(in.(mSet,correlate))
				continue
			end
			posMcounter += 1 #should be here to aovid correlated ones having a positive first than be removed
			M[mSet][:pos] = posMcounter
			
			tempmpm = []
			nowM = M[mSet][:data]
			if E[:str] == "D"
				for c in eachcol(nowM)
					push!(tempmpm,sum(c.*E[:iVarStr].*c))
				end
				M[mSet][:Mp] = map(i -> transpose(nowM[:,i].*E[:iVarStr]), axes(nowM, 2))
			else
				for c in eachcol(nowM)
					push!(tempmpm,dot(c,c))
				end
				M[mSet][:Mp] = map(i -> transpose(nowM[:,i]), axes(nowM, 2))
			end			

			M[mSet][:mpm] = tempmpm

			nowM = 0
			tempmpm = 0
			
			#summary statistics
			M[mSet][:lhs] = zeros(M[mSet][:dims][2])
			M[mSet][:rhs] = zeros(M[mSet][:dims][2])
			if mSet in keys(summaryStat)
                       		M[mSet][:lhs] .= isa(summaryStat[mSet].v,Array{Float64,1}) ? inv.(summaryStat[mSet].v) : inv.(diag(summaryStat[mSet].v))
				M[mSet][:rhs] .= isa(summaryStat[mSet].v,Array{Float64,1}) ? inv.(summaryStat[mSet].v) .* (summaryStat[mSet].m)  : inv.(diag(summaryStat[mSet].v)) .* (summaryStat[mSet].m)
				####Deal with N(0,0)
				M[mSet][:lhs][isinf.(M[mSet][:lhs])].= 0.0
				M[mSet][:rhs][isnan.(M[mSet][:rhs])].= 0.0
			end
				
			beta = push!(beta,zeros(Float64,1,M[mSet][:dims][2]))
			delta = push!(delta,ones(Int64,1,M[mSet][:dims][2]))
			
			stBWGR!(M,mSet,priorVCV,beta)
	
		elseif isa(zSet,Tuple)
			posZcounter += 1
			Z[zSet][:pos] = posZcounter
			Z[zSet][:levels] = first(getindex.(getindex.(Ref(Z),zSet),:levels))
			tempZ = hcat.(eachcol.(getindex.(getindex.(Ref(Z), zSet),:data))...)
			#same Z for all components in a single-trait model get only first column! Z[zSet][:data] = getindex.(tempZ,:,1)
			Z[zSet][:data] = tempZ
			Z[zSet][:str] = Z[zSet[1]][:str] 
			for d in zSet
                    	   	delete!(Z,d)
				delete!(modelInformation[eSet],d)
               		end
			u = push!(u,zeros(Float64,length(zSet),length(Z[zSet][:data])))
			###WEIGHTED SHOULD BE ADAAPTED HERE#################
			Z[zSet][:zpz]  = MatByMat.(tempZ)
			Z[zSet][:Zp]   = transpose.(tempZ)
			tempZ = 0
		else throw(ArgumentError("Could not understand the type of $zSet in Z"))
			
		end
		#println("KEYS OF Z: $(keys(Z))")
		#println("ZZZZZ after set: $Z")
	end

        for mSet ∈ keys(M)
		if haskey(priorVCV,mSet)
                	varBeta[mSet] = [priorVCV[mSet].v for i in 1:M[mSet][:nVarCov]]
		else
			varBeta[mSet] = [0.05 for i in 1:M[mSet][:nVarCov]]
		end
        end

	
		
	
#	for zSet in collect(keys(Z))[(!in).(keys(Z),Ref(keys(priorVCV)))]
#		posZcounter += 1
#		Z[zSet][:pos] = posZcounter
#		printstyled("No prior was provided for $pSet, but it was not included in the data. It will be made uncorrelated with default priors\n"; color = :green)		
#		tempzpz = []
#		nowZ = Z[zSet][:data]
#		for c in eachcol(nowZ)
#			push!(tempzpz,c'c)					
#		end
#		Z[zSet][:Zp]  = transpose(nowZ)						
#		Z[zSet][:zpz] = tempzpz
#		Z[zSet][:lhs] = zeros(size(nowZ,2))
#		Z[zSet][:rhs] = zeros(size(nowZ,2))
#		if zSet in keys(summaryStat)
#                	Z[zSet][:lhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) : inv.(diag(summaryStat[zSet].v))
#                        Z[zSet][:rhs] .= isa(summaryStat[zSet].v,Array{Float64,1}) ? inv.(summaryStat[zSet].v) .* (summaryStat[zSet].m)  : inv.(diag(summaryStat[zSet].v)) .* (summaryStat[zSet].m)
#                end
#	end
end


#main sampler
function getMME!(Y,X,Z,M,E,blocks,priorVCV,summaryStat,modelInformation,outPut) #maybe later use modelInformation
		
        #some info
	nRand = length(Z)
	nData = size(Y,1)
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

	println("E $E")
	println("modelInformation $modelInformation")
	#println("varE $varE")

	###################################
	
	
	######## 
	#X and b
	########

	b = []
	posXcounter = 0
	
	if isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Symbol
		println("model is a single-trait model")
		for eSet in keys(E)
			MMEX!(X,b,posXcounter,eSet,E,blocks,modelInformation,summaryStat)
		end
	elseif isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Tuple
		println("model is a multi-trait model where measurements/observations are from the same individuals")
		for eSet in keys(E)
			MMEX!(X,b,posXcounter,eSet,E,blocks,modelInformation,summaryStat)
		end
	elseif !isequal(length(collect(keys(E))),1) && all(typeof.(collect(keys(E))) .<: Symbol)
		for eSet in keys(E)
			MMEX!(X,b,posXcounter,eSet,E,blocks,modelInformation,summaryStat)
		end
		println("model is a multi-population model where measurements/observations are from different individuals")
	else 	throw(ArgumentError("Could not understand the type of your model"))

	end

	

	###Allow no fixed effects
	#isempty(keys(X)) ? b = [] : b = zeros(sum(getindex.(getindex.(Ref(X), keys(X)),:nCol)))
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)

	###################################
	
	### 
	#Z and u
	###
	
	u = []	
	posZcounter = 0

	
	if isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Symbol
		println("model is a single-trait model")
		for eSet in keys(E)
			MMEZ!(Z,u,varU,posZcounter,eSet,E,priorVCV,modelInformation,summaryStat)
		end
		varCovZ!(Z,priorVCV)
	elseif isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Tuple
		println("model is a multi-trait model where measurements/observations are from the same individuals")
#		for eSet in keys(E)
#			MMEZ!(Z,u,varU,posZcounter,eSet,priorVCV,modelInformation,summaryStat)
#		end
	elseif !isequal(length(collect(keys(E))),1) && all(typeof.(collect(keys(E))) .<: Symbol)
		for eSet in keys(E)
			MMEZ!(Z,u,varU,posZcounter,eSet,priorVCV,modelInformation,summaryStat)
		end
		println("model is a multi-population model where measurements/observations are from different individuals")
	else 	throw(ArgumentError("Could not understand the type of your model"))

	end

																		

	#ADD MARKERS
	# read map file and make regions
																		
	beta = []
	delta = []
	posMcounter = 0

	if isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Symbol
		println("model is a single-trait model")
		for eSet in keys(E)
			MMEM!(M,beta,varBeta,posMcounter,eSet,E,priorVCV,modelInformation,summaryStat)
		end
		varCovM!(M,priorVCV)
	elseif isequal(length(collect(keys(E))),1) && typeof(collect(keys(E))[]) <: Tuple
		println("model is a multi-trait model where measurements/observations are from the same individuals")
#		for eSet in keys(E)
#			MMEM!(M,beta,varBeta,posMcounter,eSet,E,priorVCV,modelInformation,summaryStat)
#		end
	elseif !isequal(length(collect(keys(E))),1) && all(typeof.(collect(keys(E))) .<: Symbol)
		for eSet in keys(E)
#			MMEM!(M,beta,varBeta,posMcounter,eSet,E,priorVCV,modelInformation,summaryStat)
		end
		println("model is a multi-population model where measurements/observations are from different individuals")
	else 	throw(ArgumentError("Could not understand the type of your model"))
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

	[inOut.outMCMC(outPut,"b_$eSet",levelsX) for eSet in keys(E) if isa(eSet,Symbol)]
	[inOut.outMCMC(outPut,"b_$eSet",hcat(["$(l)_".*String.(hcat(eSet...)) for l in levelsX]...)) for eSet in keys(E) if isa(eSet,Tuple{Vararg{Symbol}})]
	
	
	#check for correlated RE
 #       for zSet in keys(Z)
#		if isa(zSet, Symbol)
#			nameRE_VCV = String(zSet)
#			inOut.outMCMC(outPut,"u$zSet",[Z[zSet][:levels]])
#			inOut.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
#		elseif isa(zSet, Expr)
#			nameRE_VCV = join(zSet.args)[2:end]
#			inOut.outMCMC(outPut,"u$zSet",[Z[zSet][:levels]])
#			inOut.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
#		elseif isa(zSet, Tuple)
#			nameRE_VCV =  join(String.(vcat(zSet...)),"_").*hcat(["_$i" for i in 1:(length(zSet)^2)]...)
#			for z in zSet
 #  				inOut.outMCMC(outPut,"u$z",[Z[zSet][:levels]])
#			end
#			inOut.outMCMC(outPut,"varU$zSet",[nameRE_VCV]) #[] to have it as one row
#		end
#	end	

	[[inOut.outMCMC(outPut,"u$zSet$eSet",[Z[zSet][:levels]]) for zSet in keys(Z) if (isa(zSet,Union{Symbol,Expr}) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single trait	
	
	#this prints all u into one file ID_in1,Dam_ind1,ID_ind2,Dam_ind2.... order
	[[inOut.outMCMC(outPut,"u"*join([String.(zSet)...])*"$eSet",hcat(vcat([Array([String.(zSet)...].*i) for i in Z[zSet][:levels]]...)...)) for zSet in keys(Z) if (isa(zSet,Tuple) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Symbol)] #single-trait multiple comp
	#this prints all u into seperate files
	#[[[inOut.outMCMC(outPut,"u$zSet$eSet",[Z[zSet][:levels]]) for z in zSet] for zSet in keys(Z) if (isa(zSet,Tuple) && in(zSet,keys(modelInformation[eSet])))] for eSet in keys(E) if isa(eSet,Tuple)] #single-trait multiple comp

	
	#[(inOut.outMCMC(outPut,"varU$zSet",[String(zSet)]) for zSet in keys(Z) if isa(zSet,Symbol)) for eSet in keys(E) if isa(eSet,Symbol)] #[] to have it as one row
	#[(inOut.outMCMC(outPut,"varU$zSet",[join(zSet.args)[2:end]]) for zSet in keys(Z) if isa(zSet,Expr)) for eSet in keys(E) if isa(eSet,Symbol)] #[] to have it as one row	
	#[(inOut.outMCMC(outPut,"varU$zSet",[join(String.(vcat(zSet...)),"_").*hcat(["_$i" for i in 1:(length(zSet)^2)]...)])) for zSet in keys(Z) for eSet in keys(E) if isa(eSet,Tuple)]
	
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
	

	
	[inOut.outMCMC(outPut,"varE_$eSet",["e_$eSet"]) for eSet in keys(E) if isa(eSet,Symbol)]
	[inOut.outMCMC(outPut,"varE_$eSet",hcat(["$(c[1])_$(c[2])" for c in collect(Iterators.product(1:length(eSet), 1:length(eSet)))]...)) for eSet in keys(E) if isa(eSet,Tuple{Vararg{Symbol}})]

	
	##########
	
	X  = myUnzip(X)
	Z  = myUnzip(Z)
	M  = myUnzip(M)
	E  = myUnzip(E) #(;E...)
	
	return modelInformation,ycorr, nData, E, varE, X, b, Z, u, varU, M,  beta, varBeta, delta
	
end

