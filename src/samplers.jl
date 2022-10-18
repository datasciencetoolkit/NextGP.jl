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

export runSampler

#main sampler
function runSampler(iA,Y,X,Z,levelDict,blocks,chainLength,burnIn,outputFreq,priorVCV,M,paths2maps,rS,outPut)
		
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	nRand = length(Z)
	nColEachZ    = OrderedDict(z => size(Z[z],2) for z in keys(Z))
	nData = length(Y)
	nMarkerSets = length(M)
	nMarkers    = OrderedDict(m => size(M[m],2) for m in keys(M))

        #initial computations and settings
	ycorr = deepcopy(Y)
	
	### X and b
	levelsX = levelDict[:levelsFE]
	
	#==BLOCK FIXED EFFECTS.
	Order of blocks is as definde by the user
	Order of variables within blocks is always the same as in the model definition, not defined by the user in each block.
	==#
	for b in blocks
		getThese = intersect(collect(keys(X)), b)
		X[Tuple(getThese)] = hcat(getindex.(Ref(X), getThese)...)
		levelsX[Tuple(getThese)] = vcat(getindex.(Ref(levelsX), getThese)...)
		for d in getThese
			delete!(X,d)
			delete!(levelsX,d)
		end
	end
	
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	
	#not a dictionary anymore, and consistent with possible new order.
	levelsX = hcat(vcat([isa(value,String) ? value : vcat(value...) for (key, value) in levelsX]...)...)
	
	#Key positions of variablese and blocks for speed. b is an array of arrays.
        XKeyPos = OrderedDict{Any,Int64}()
        [XKeyPos[collect(keys(X))[i]]=i for i in 1:length(keys(X))]
	        
	##make iXpX, Z', zpz (for uncor)
        iXpX = deepcopy(X)
        for x in keys(X)
		XpX = X[x]'X[x]
		if isa(XpX,Matrix{Float64}) 
			XpX += Matrix(I*minimum(abs.(diag(XpX)./10000)),size(XpX))
			#XpX += Matrix(I*minimum(abs.(diag(XpX)./size(X[x],1))),size(XpX))
		end
               	iXpX[x] = inv(XpX)
        end
	
        ##make b and u arrays
        b = Array{Array{Float64, 1},1}(undef,0)
        ##counts columns per effect
        nColEachX = []
        for xSet in keys(X)
                nCol = size(X[xSet],2)
                push!(b,fill(0.0,nCol))
                nColEachX = push!(nColEachX,nCol)
        end

	#set up for E.
						
	#no inverse implemented yet!
	if haskey(priorVCV,:e)	
		if isempty(priorVCV[:e][1]) || priorVCV[:e][1]=="I" 
				printstyled("prior var-cov structure for \"e\" is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				strE = Matrix(1.0I,nData,nData)
				priorVCV[:e] = ("I",priorVCV[:e][2])
		elseif priorVCV[:e][1]=="D"
				strE = D ##no inverse  yet
				error("var-cov structure \"D\" has not been implemented yet")
				printstyled("prior var-cov structure for \"e\" is \"D\". User provided \"D\" matrix (d_ii = 1/w_ii) will be used\n"; color = :green)
		else 
				error("provide a valid prior var-cov structure (\"I\", \"D\" or leave it empty \"[]\") for \"e\" ")
		end
	else	
		printstyled("prior var-cov for \"e\" is fully  empty. An identity matrix will be used with an arbitrary variance of 100\n"; color = :green)
		strE = Matrix(1.0I,nData,nData)
		varE_prior = 100
		#just add to priors
		priorVCV[:e] = ("I",varE_prior)
	end
								
	#parameters for priors
        dfE = 4.0
	dfDefault = 4.0
 
	       
	if priorVCV[:e][2]==0.0
		priorVCV[:e][2]  = 0.0005
       		scaleE     = 0.0005
        else
       		scaleE    = priorVCV[:e][2]*(dfE-2.0)/dfE    
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
			if isempty(priorVCV[zSet][1]) || priorVCV[zSet][1]=="I" 
				printstyled("prior var-cov structure for $zSet is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				iVarStr[zSet] = Matrix(1.0I,nCol,nCol)
			elseif priorVCV[zSet][1]=="A"
				iVarStr[zSet] = iA
				printstyled("prior var-cov structure for $zSet is A. Computed A matrix (from pedigree file) will be used\n"; color = :green)
			else 	iVarStr[zSet] = inv(priorVCV[zSet][1])
			end
			varU_prior[zSet] = priorVCV[zSet][2]
		else	
			printstyled("prior var-cov for $zSet is empty. An identity matrix will be used with an arbitrary variance of 100\n"; color = :green)
		iVarStr[zSet] = Matrix(1.0I,nCol,nCol)
		varU_prior[zSet] = 100
		priorVCV[zSet] = ("I",varU_prior[zSet])
		end
        end

	#df, shape, scale...															
	
	dfZ = Dict{Any,Any}()	
	for zSet ∈ keys(zpz)
		dfZ[zSet] = 3.0+size(priorVCV[zSet][2],1)
	end
																
	scaleZ = Dict{Any,Any}()
        for zSet in keys(zpz)
                nZComp = size(priorVCV[zSet][2],1)
		#priorVCV[zSet][2] is a temporary solution
		nZComp > 1 ? scaleZ[zSet] = priorVCV[zSet][2].*(dfZ[zSet]-nZComp-1.0)  : scaleZ[zSet] = priorVCV[zSet][2]*(dfZ[zSet]-2.0)/dfZ[zSet] #I make float and array of float														
        end

												
        ####
																					

	#ADD MARKERS
	# read map file and make regions
																		
	############priorVCV cannot be empty for markers, currently!!																	

	#key positions for each effect in beta, for speed. Order of matrices in M are preserved here.
        BetaKeyPos = OrderedDict{Any,Any}()
        for mSet in keys(M)
                pos = findall(mSet.==collect(keys(M)))[]
                BetaKeyPos[mSet] = pos
        end


	#make mpm

	Mp = OrderedDict{Any,Any}()
       	mpm = OrderedDict{Any,Any}()
		
	regionArray = OrderedDict{Any,Array{UnitRange{Int64},1}}()	
	
	for pSet ∈ keys(filter(p -> p.first!=:e, priorVCV)) # excluding :e keys(priorVCV)
		corEffects = []
		corPositions = []
		#symbol :M1 or expression
		if isa(pSet,Symbol) && in(pSet,keys(M))
			tempmpm = []
			nowM = M[pSet]
			for c in eachcol(nowM)
				push!(tempmpm,BLAS.dot(c,c))
			end
			mpm[pSet] = tempmpm
			theseRegions = prep2RegionData(outPut,pSet,paths2maps[pSet],rS[pSet])
		        regionArray[pSet] = theseRegions
		#tuple of symbols (:M1,:M2)
		elseif (isa(pSet,Tuple{Vararg{Symbol}})) && all((in).(pSet,Ref(keys(M)))) #if all elements are available # all([pSet .in Ref(keys(M))])
			correlate = collect(pSet)
			for pSubSet in correlate
				push!(corEffects,pSubSet)
				push!(corPositions,findall(pSubSet.==keys(M))[])
			end
			if issubset(corEffects,collect(keys(M)))
				tempM = hcat.(eachcol.(getindex.(Ref(M), (pSet)))...)
				rS[pSet] = rS[first(pSet)]
				for d in corEffects
                       			delete!(M,d)
					delete!(BetaKeyPos,d)
					delete!(rS,d)
               			end
				BetaKeyPos[pSet] = corPositions
				M[pSet]   = tempM
				mpm[pSet] = MatByMat.(tempM)
				Mp[pSet]  = transpose.(tempM)
				tempM = 0
				nowMap = first(pSet)		#should throw out error if sets have different lengths! implement it here!
				theseRegions = prep2RegionData(outPut,pSet,paths2maps[nowMap],rS[pSet]) ###first data
                		regionArray[pSet] = theseRegions
			end
		end
	end
	
	for pSet in collect(keys(M))[(!in).(keys(M),Ref(keys(priorVCV)))]
		printstyled("No prior was provided for $pSet, but it was included in the data. It will be made uncorrelated with default priors\n"; color = :green)		
		tempmpm = []
		nowM = M[pSet]
		for c in eachcol(nowM)
			push!(tempmpm,BLAS.dot(c,c))
		end
		mpm[pSet] = tempmpm
		theseRegions = prep2RegionData(outPut,pSet,paths2maps[pSet],rS[pSet])
		regionArray[pSet] = theseRegions
	end
	
	#pos for individual marker set
	BetaKeyPos4Print = OrderedDict(vcat([isa(k,Symbol) ? k => v : collect(k) .=> collect(v) for (k,v) in BetaKeyPos]...))

	nRegions  = OrderedDict(mSet => length(regionArray[mSet]) for mSet in keys(regionArray))
	
	dfM = Dict{Any,Any}()	
	for mSet ∈ keys(mpm)
		dfM[mSet] = 3.0+size(priorVCV[mSet],1)
	end


	scaleM = Dict{Any,Any}()
        for mSet in keys(mpm)
                nMComp = size(priorVCV[mSet],1)
                nMComp > 1 ? scaleM[mSet] = priorVCV[mSet].*(dfM[mSet]-nMComp-1.0)  : scaleM[mSet] = priorVCV[mSet]*(dfM[mSet]-2.0)/dfM[mSet] #I make float and array of float
        end
	
	
	#storage
	u = zeros(Float64,nRand,maximum(vcat([0,collect(values(nColEachZ))]...))) #zero is for max to work when no random effect is present #can allow unequal length! Remove tail zeros for printing....

	varU = deepcopy(varU_prior) #for storage

	beta = zeros(Float64,nMarkerSets,maximum(vcat([0,collect(values(nMarkers))]...))) #zero is for max to work when no SNP data is present #can allow unequal length! Remove tail zeros for printing....

        varBeta = OrderedDict{Any,Any}()
        for mSet in keys(mpm)
                varBeta[mSet] = [priorVCV[mSet] for i in 1:length(regionArray[mSet])] #later, direct reference to key when varM_prior is a dictionary
        end

	#summarize analysis
	summarize = DataFrame(Effect=Any[],Type=Any[],Str=Any[],df=Any[],scale=Any[])
	
	for zSet in keys(zpz)
		if zSet ∈ keys(priorVCV)
			str = priorVCV[zSet][1]
			#value = priorVCV[zSet][2]
		else 
			str = "I"
		     	#value = varU_prior[zSet]
		end
	push!(summarize,[zSet,"Random",str,dfZ[zSet],scaleZ[zSet]])
	end
	for mSet in keys(mpm)
		if mSet ∈ keys(priorVCV)
			str = "$(nRegions[mSet]) block(s)"
			#value = priorVCV[mSet]
		else #### later, handel this above, when dealing with priorVCV is allowed to be empty
			str = "WG(I)"
		     	#value = 0.001
		end
	push!(summarize,[mSet,"Random (Marker)",str,dfM[mSet],scaleM[mSet]])
	end
	push!(summarize,["e","Random",priorVCV[:e][1],dfE,scaleE])						

	println("\n ---------------- Summary of analysis ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,nosubheader=true,alignment=:l)


	#########make MCMC output files.
	IO.outMCMC(outPut,"b",levelsX)
	
	#check for correlated RE
        for i in 1:length(levelDict[:levelsRE])
		levRE = hcat(vcat(collect(values(levelDict[:levelsRE]))[i]...)...)
		IO.outMCMC(outPut,"u$i",levRE)
		isa(collect(keys(levelDict[:levelsRE]))[i], Symbol) ? nameRE_VCV = String(collect(keys(levelDict[:levelsRE]))[i]) : nameRE_VCV = join(collect(keys(levelDict[:levelsRE]))[i].args)[2:end]
#		IO.outMCMC(outPut,"varU$i",[join(collect(keys(levelDict[:levelsRE]))[i],"_")])
		IO.outMCMC(outPut,"varU$i",[nameRE_VCV]) #[] to have it as one row
	end	
	
	#arbitrary marker names
	for mSet in keys(BetaKeyPos4Print)
   		IO.outMCMC(outPut,"beta$mSet",hcat(["rs_$i" for i in 1:nMarkers[mSet]]...))
		isa(mSet, Symbol) ? nameM_VCV = ["reg_$r" for r in 1:nRegions[mSet]] : nameRE_VCV = vcat([["reg_$i" for j in 1:length(mSet)^2] for i in 1:nRegions[mSet]]...)
		IO.outMCMC(outPut,"var$mSet",[nameM_VCV]) #[] to have it as one row
        end
	

	IO.outMCMC(outPut,"varE",["e"])
	##########


	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
	sleep(0.1)
	
		#sample residual variance
	       	varE = sampleVarE(dfE,scaleE,ycorr,nData)
		
		#sample fixed effects
	        sampleX!(X,b,iXpX,nFix,nColEachX,XKeyPos,ycorr,varE)
	
		#sample random effects
	        sampleZandZVar!(iVarStr,Z,Zp,u,zpz,uKeyPos,nColEachZ,ycorr,varE,varU,scaleZ,dfZ)	

		#sample marker effects and variances
	        sampleMandMVar_view!(M,Mp,beta,mpm,nMarkerSets,BetaKeyPos,regionArray,nRegions,ycorr,varE,varBeta,scaleM,dfM)
               		
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
                                IO.outMCMC(outPut,"beta$mSet",beta[BetaKeyPos4Print[mSet],:]')
                        end
			for pSet in keys(mpm)
				IO.outMCMC(outPut,"var".*String(pSet),varBeta[pSet]')
			end
		end
	end
end


#Sampling fixed effects

function sampleX!(X,b,iXpX,nFix,nColEachX,keyX,ycorr,varE)
        #block for each effect
        for xSet in keys(X)
		pos = keyX[xSet]	
                if nColEachX[pos] == 1
			ycorr    .+= X[xSet].*b[pos]
			rhs      = X[xSet]'*ycorr
			meanMu   = iXpX[xSet]*rhs			
                        b[pos] .= rand(Normal(meanMu[],sqrt((iXpX[xSet]*varE))[]))
			ycorr    .-= X[xSet].*b[pos]
                else	
			ycorr    .+= X[xSet]*b[pos]
                        rhs      = X[xSet]'*ycorr
                        meanMu   = iXpX[xSet]*rhs
			b[pos] .= rand(MvNormal(vec(meanMu),convert(Array,Symmetric(iXpX[xSet]*varE))))
			ycorr    .-= X[xSet]*b[pos]
                end
        end
end


function sampleU(iMat,pos,ZComp,ZpComp,zpzComp,varE,varUComp,uVector,ycorr)
	uVec = deepcopy(uVector)
	λz = varE/varUComp
	Yi = ZpComp*ycorr #computation of Z'ycorr for ALL  rhsU
	nCol = length(zpzComp)
	for i in 1:nCol
        	uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
		rhsU = Yi[i] - λz*dot(iMat[:,i],uVec)
                lhsU = zpzComp[i] + (iMat[i,i]*λz)[1]
		invLhsU = 1.0/lhsU
                meanU = invLhsU*rhsU
                uVec[i] = rand(Normal(meanU,sqrt(invLhsU*varE)))
        end
	return uVec
end


function sampleZandZVar!(iStrMat,ZMat,ZpMat,u,zpzMat,keyU,nCols,ycorr,varE,varU,scaleZ,dfZ)
        #for each random effect
        for zSet in keys(zpzMat)
		if isa(zSet,Tuple)
			uPos = keyU[zSet]
			nowZp = ZpMat[zSet] ###
			error("correlated random effects are not allowed")
		elseif isa(zSet,Symbol) || isa(zSet,Expr)
                	uPos = keyU[zSet]
			ycorr .+= ZMat[zSet]*u[uPos,1:nCols[zSet]]
                	u[uPos,1:nCols[zSet]]  .= sampleU(iStrMat[zSet],uPos,ZMat[zSet],ZpMat[zSet],zpzMat[zSet],varE,varU[zSet],u[uPos,1:nCols[zSet]],ycorr)
			ycorr .-= ZMat[zSet]*u[uPos,1:nCols[zSet]]		
			varU[zSet] = sampleVarU(iStrMat[zSet],scaleZ[zSet],dfZ[zSet],u[uPos,1:nCols[zSet]])
       		end
		
	 end
end

function sampleMandMVar_view!(MMat,MpMat,beta,mpmMat,nMSet,keyBeta,regionsMat,regions,ycorr,varE,varBeta,scaleM,dfM)
        #for each marker set
        for mSet in keys(mpmMat)
		if isa(mSet,Tuple)
			betaPos = keyBeta[mSet]
			nowMp = MpMat[mSet] ###
			for r in 1:regions[mSet]
				theseLoci = regionsMat[mSet][r]
				regionSize = length(theseLoci)
				invB = inv(varBeta[mSet][r])
				for locus in theseLoci
					RHS = zeros(size(invB,1))	
					ycorr .+= MMat[mSet][locus]*beta[betaPos,locus]					
					RHS = (nowMp[locus]*ycorr)./varE
					invLHS::Array{Float64,2} = inv((mpmMat[mSet][locus]./varE) .+ invB)
					meanBETA::Array{Float64,1} = invLHS*RHS
					beta[betaPos,locus] = rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS))))
					ycorr .-= MMat[mSet][locus]*beta[betaPos,locus]	
				end
				varBeta[mSet][r] = sampleVarCovBeta(scaleM[mSet],dfM[mSet],beta[betaPos,theseLoci],regionSize)
			end	
		else
		println("-------------------")
			nowM = MMat[mSet]
                	betaPos = keyBeta[mSet]
			local rhs::Float64
			local lhs::Float64
			local meanBeta::Float64
                	for r in 1:regions[mSet]
                        	theseLoci = regionsMat[mSet][r]
                        	regionSize = length(theseLoci)
                        	lambda = varE/(varBeta[mSet][r])
                        	for locus in theseLoci
					BLAS.axpy!(beta[betaPos,locus],nowM[:,locus],ycorr)
                                	rhs = BLAS.dot(nowM[:,locus],ycorr)
	                               	lhs = mpmMat[mSet][locus] + lambda
                                	meanBeta = lhs\rhs
                                	beta[betaPos,locus] = sampleBeta(meanBeta, lhs, varE)
                                	BLAS.axpy!(-1.0*beta[betaPos,locus],nowM[:,locus],ycorr)
                       		end
                        	varBeta[mSet][r] = sampleVarBeta(scaleM[mSet],dfM[mSet],beta[betaPos,theseLoci],regionSize)
                	end
       		end
	 end
end


#sample marker effects
function sampleBeta(meanBeta, lhs, varE)
    return rand(Normal(meanBeta,sqrt(lhs\varE)))
end

#sample random effects' variances (new U)
function sampleVarU(iMat,scale_ranVar,df_ranVar,effVec)
	n = size(iMat,2)
	return (scale_ranVar*df_ranVar + effVec'*iMat*effVec)/rand(Chisq(df_ranVar + n))
end

#sample marker variances
function sampleVarBeta(scalem,dfm,whichLoci,regionSize)
	return (scalem*dfm + BLAS.dot(whichLoci,whichLoci)) / rand(Chisq(dfm + regionSize))
end

function sampleVarCovBeta(scalem,dfm,whichLoci,regionSize)
	Sb = whichLoci*whichLoci'
	return rand(InverseWishart(dfm + regionSize, scalem + Sb))
end

#Sample residual variance
function sampleVarE(df_e,S_e,yCorVec,nRecords)
    return (df_e*S_e + BLAS.dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords))
end



end
