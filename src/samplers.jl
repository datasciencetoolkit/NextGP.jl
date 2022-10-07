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
function runSampler(iA,Y,X,Z,levelDict,chainLength,burnIn,outputFreq,priorVCV,M,paths2maps,rS,outPut)
	
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	nRand = length(Z)
	nColEachZ    = OrderedDict(z => size(Z[z],2) for z in keys(Z))
	nData = length(Y)
	nMarkerSets = length(M)
	nMarkers    = [size(M[m],2) for m in keys(M)]

        #initial computations and settings
	ycorr = deepcopy(Y)
	
	### X and b
	
	#Key positions for speed. Old positions, before blocking!
        XKeyPos = OrderedDict{Any,Int64}()
        [XKeyPos[collect(keys(X))[i]]=i for i in 1:length(keys(X))]

	
	levelsX = levelDict[:levelsFE]
	
	#BLOCK FIXED EFFECTS
	for b in blocks
		getThese = intersect(collect(keys(X)), b)
		X[Tuple(getThese)] = hcat(getindex.(Ref(X), getThese)...)
		levelsX[Tuple(getThese)] = vcat(getindex.(Ref(levelsX), getThese)...)
		for d in getThese
			delete!(X,d)
			delete!(levelsX,d)
		end
	end
	
	#not a dictionary anymore, and consistent with possible new order.
	levelsX = hcat(vcat([isa(value,String) ? value : vcat(value...) for (key, value) in levelsX]...)...)
	        
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
#                println(xSet)
                nCol = size(X[xSet],2)
                push!(b,fill(0.0,nCol))
                nColEachX = push!(nColEachX,nCol)
        end

	#set up for E
#	isempty(priorVCV["e"][1]) ? strE = Matrix(1.0I,nData,nData) : strE = priorVCV["e"][1]
#	isempty(priorVCV["e"][2]) ? varE_prior = 100 : varE_prior = priorVCV["e"][2]
	
	#no inverse implemented yet!
	if haskey(priorVCV,"e")	
		if isempty(priorVCV["e"][1]) || priorVCV["e"][1]=="I" 
				printstyled("prior var-cov structure for \"e\" is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				strE = Matrix(1.0I,nData,nData)
				priorVCV["e"] = ("I",priorVCV["e"][2])
		elseif priorVCV["e"][1]=="D"
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
		priorVCV["e"] = ("I",varE_prior)
	end
				

	#parameters for priors
        dfE = 4.0
	dfDefault = 4.0
 
	       
	if priorVCV["e"][2]==0.0
		priorVCV["e"][2]  = 0.0005
       		scaleE     = 0.0005
        else
       		scaleE    = priorVCV["e"][2]*(dfE-2.0)/dfE    
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
       	zpz = OrderedDict{Any,Any}()
		
	corZ = OrderedDict{Any,Any}()
	corZPos = OrderedDict{Any,Any}()
											
	for pSet ∈ keys(priorVCV)
		corEffects = []
		corPositions = []
		if typeof(pSet)==Tuple{String, String}
			if pSet ∈ keys(Z)
				tempzpz = []
				nowZ = Z[pSet]
				for c in eachcol(nowZ)
					push!(tempzpz,c'c)					
					# push!(tempzpz,BLAS.dot(c,c))
				end
				Zp[pSet]  = transpose(Z[pSet])						
				zpz[pSet] = tempzpz
			end
		elseif  issubset(pSet,keys(Z))
			correlate = collect(pSet)
			for pSubSet in correlate
				push!(corEffects,pSubSet)
				push!(corPositions,findall(pSubSet.==keys(Z))[])
			end
			if issubset(corEffects,collect(keys(Z)))
				corZPos[pSet] = corPositions
				corZ[pSet] = corEffects
				tempZ = hcat.(eachcol.(getindex.(Ref(Z), (pSet)))...)
				for d in corEffects
                       			delete!(Z,d)
               			end
				Z[pSet]   = tempZ
				zpz[pSet] = MatByMat.(tempZ)
				Zp[pSet]  = transpose.(tempZ)
				tempZ = 0
			end
		end
	end
	
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
        BetaKeyPos = OrderedDict{String,Int64}()
        for mSet in keys(M)
                pos = findall(mSet.==collect(keys(M)))[]
                BetaKeyPos[mSet] = pos
        end


	#make mpm

	Mp = OrderedDict{Any,Any}()
       	mpm = OrderedDict{Any,Any}()
		
	corM = OrderedDict{Any,Any}()
	corMPos = OrderedDict{Any,Any}()
	regionArray = OrderedDict{Any,Array{UnitRange{Int64},1}}()	

	for pSet ∈ keys(priorVCV)
		corEffects = []
		corPositions = []
		if typeof(pSet)==String
			println("$pSet is univariate")
			if pSet ∈ keys(M)
				println("univariate mpm for $pSet")
				tempmpm = []
				nowM = M[pSet]
				for c in eachcol(nowM)
					push!(tempmpm,BLAS.dot(c,c))
				end
				mpm[pSet] = tempmpm
				theseRegions = prep2RegionData(outPut,pSet,paths2maps[pSet],rS[pSet])
		                regionArray[pSet] = theseRegions
			end
		elseif  issubset(pSet,keys(M))
			println("$pSet will be correlated")
			correlate = collect(pSet)
			for pSubSet in correlate
				println(pSubSet)
				push!(corEffects,pSubSet)
				push!(corPositions,findall(pSubSet.==keys(M))[])
			end
			if issubset(corEffects,collect(keys(M)))
				corMPos[pSet] = corPositions
				corM[pSet] = corEffects
				tempM = hcat.(eachcol.(getindex.(Ref(M), (pSet)))...)
				for d in corEffects
                       			delete!(M,d)
               			end
				M[pSet]   = tempM
				mpm[pSet] = MatByMat.(tempM)
				Mp[pSet]  = transpose.(tempM)
				tempM = 0
				nowMap = first(pSet)					 #should throw out error if sets have different lengths! implement it here!
				theseRegions = prep2RegionData(outPut,pSet,paths2maps[nowMap],rS[nowMap]) ###first data
                		regionArray[pSet] = theseRegions
			end
		end
	end  	

	nRegions  = OrderedDict(mSet => length(regionArray[mSet]) for mSet in keys(regionArray))
	println("number of regions: $(values(nRegions))")


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

	beta = zeros(Float64,nMarkerSets,maximum(vcat([0,nMarkers]...))) #zero is for max to work when no SNP data is present #can allow unequal length! Remove tail zeros for printing....
#	vcovBeta = fill(Matrix(Diagonal(varM)),maximum(nRegions)) #can allow unequal length! Remove tail zeros for printing....

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
	push!(summarize,[zSet,"Z",str,dfZ[zSet],scaleZ[zSet]])
	end
	for mSet in keys(mpm)
		if mSet ∈ keys(priorVCV)
			str = nRegions[mSet]
			#value = priorVCV[mSet]
		else #### later, handel this above, when dealing with priorVCV is allowed to be empty
			str = "WG(I)"
		     	#value = 0.001
		end
	push!(summarize,[mSet,"M",str,dfM[mSet],scaleM[mSet]])
	end
	push!(summarize,["e","Res",priorVCV["e"][1],dfE,scaleE])
	println("\n ---------------- Summary of analysis ---------------- \n")
	pretty_table(summarize, tf = tf_markdown, show_row_number = false,nosubheader=true,alignment=:l)


	#########make MCMC output files.
	IO.outMCMC(outPut,"b",levelsX)

	#check for correlated RE
        for i in 1:length(levelDict[:levelsRE])
		nameRE = hcat(vcat(collect(values(levelDict[:levelsRE]))[i]...)...)
		IO.outMCMC(outPut,"u$i",nameRE)
		IO.outMCMC(outPut,"varU$i",[join(collect(keys(levelDict[:levelsRE]))[i],"_")])
	end	
	

	IO.outMCMC(outPut,"varE",["varE"])
	##########


	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
	sleep(0.1)
	
		#sample residual variance
	       	varE = sampleVarE(dfE,scaleE,ycorr,nData)
		
		#sample fixed effects
	        sampleX!(X,b,iXpX,nFix,nColEachX,XKeyPos,ycorr,varE)
	
		#sample random effects
	        sampleZandZVar!(iVarStr,Z,Zp,corZ,corZPos,u,zpz,uKeyPos,nColEachZ,ycorr,varE,varU,scaleZ,dfZ)	

		#sample marker effects and variances
	        sampleMandMVar_view!(M,Mp,corM,corMPos,beta,mpm,nMarkerSets,BetaKeyPos,regionArray,nRegions,ycorr,varE,varBeta,scaleM,dfM)
               		
        	#print
		if iter in these2Keep
			IO.outMCMC(outPut,"b",vcat(b...)') ### currently no path is provided!!!!
			IO.outMCMC(outPut,"varE",varE)
			
			for zSet in keys(uKeyPos)
                                IO.outMCMC(outPut,"u$(uKeyPos[zSet])",u[uKeyPos[zSet],1:nColEachZ[zSet]]')
                        end
			for pSet in keys(zpz)
				IO.outMCMC(outPut,"varU$(uKeyPos[pSet])",varU[pSet]) #join values for multivariate in uKeyPos[pSet])
			end

			for mSet in keys(BetaKeyPos)
                                IO.outMCMC(outPut,"beta$mSet",beta[BetaKeyPos[mSet],:]')
                        end
			for pSet in keys(mpm)
				IO.outMCMC(outPut,"var".*pSet,varBeta[pSet])
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


function sampleZandZVar!(iStrMat,ZMat,ZpMat,correlatedZ,keyCorZ,u,zpzMat,keyU,nCols,ycorr,varE,varU,scaleZ,dfZ)
        #for each random effect
        for zSet in keys(zpzMat)
		if zSet in keys(correlatedZ)
			uPos = keyCorZ[zSet]
			nowZp = ZpMat[zSet] ###
			error("correlated random effects are not allowed")
		else
                	uPos = keyU[zSet]
			ycorr .+= ZMat[zSet]*u[uPos,1:nCols[zSet]]
                	u[uPos,1:nCols[zSet]]  .= sampleU(iStrMat[zSet],uPos,ZMat[zSet],ZpMat[zSet],zpzMat[zSet],varE,varU[zSet],u[uPos,1:nCols[zSet]],ycorr)
			ycorr .-= ZMat[zSet]*u[uPos,1:nCols[zSet]]		
			varU[zSet] = sampleVarU(iStrMat[zSet],scaleZ[zSet],dfZ[zSet],u[uPos,1:nCols[zSet]])
       		end
		
	 end
end

function sampleMandMVar_view!(MMat,MpMat,correlatedM,keyCorM,beta,mpmMat,nMSet,keyBeta,regionsMat,regions,ycorr,varE,varBeta,scaleM,dfM)
        #for each marker set
        for mSet in keys(mpmMat)
		if mSet in keys(correlatedM)
			betaPos = keyCorM[mSet]
			nowMp = MpMat[mSet] ###
			for r in 1:regions[mSet]
				theseLoci = regionsMat[mSet][r]
				regionSize = length(theseLoci)
				invB = inv(varBeta[mSet][r])
				for locus in theseLoci
					RHS = zeros(size(invB,1))###
				#	for m in correlatedM[mSet] 
				#		BLAS.axpy!(beta[keyBeta[m],locus],MMat[m][:,locus],ycorr) #beta pos is different than pos
				#	end
					
					ycorr .+= MMat[mSet][locus]*beta[betaPos,locus]					

					RHS = (nowMp[locus]*ycorr)./varE ### FASTEST
				#	RHS = zeros(size(invB,1)) ### for mul!
				#	mul!(RHS,nowMp[locus],ycorr./varE) ### LESS MEMORY ALLOCATION
				#	RHS = [BLAS.dot(MMat[m][:,locus],ycorr)/varE for m in correlatedM[mSet]]
					invLHS::Array{Float64,2} = inv((mpmMat[mSet][locus]./varE) .+ invB)
					meanBeta::Array{Float64,1} = invLHS*RHS
					beta[betaPos,locus] = rand(MvNormal(meanBeta,convert(Array,Symmetric(invLHS))))
				#	for m in correlatedM[mSet]
                                #                BLAS.axpy!(-1.0*beta[keyBeta[m],locus],MMat[m][:,locus],ycorr)
                                #       end
					ycorr .-= MMat[mSet][locus]*beta[betaPos,locus]	
				end
				varBeta[mSet][r] = sampleVarCovBeta(scaleM[mSet],dfM[mSet],beta[betaPos,theseLoci],regionSize)
			end	
		else
			nowM = MMat[mSet]
                	betaPos = keyBeta[mSet]
                	for r in 1:regions[mSet]
                        	theseLoci = regionsMat[mSet][r]
                        	regionSize = length(theseLoci)
                        	lambda = varE/(varBeta[mSet][r])
                        	for locus in theseLoci
                                	BLAS.axpy!(beta[betaPos,locus],view(nowM,:,locus),ycorr)
                                	rhs::Float64 = BLAS.dot(view(nowM,:,locus),ycorr)
                               		lhs::Float64 = mpmMat[mSet][locus] + lambda
                                	meanBeta::Float64 = lhs\rhs
                                	beta[betaPos,locus] = sampleBeta(meanBeta, lhs, varE)
                                	BLAS.axpy!(-1.0*beta[betaPos,locus],view(nowM,:,locus),ycorr)
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

function sampleMarkerVar!(beta,varBeta,nMSet,keyM,regions,regionsMat,scaleM,dfM)
        #for each marker set
        for mSet in keys(varBeta)
		pos = keyM[mSet]
                for r in 1:regions[pos] #dont have to compute 1000000 times, take it out
                        theseLoci = regionsMat[pos][r]
                        regionSize = length(theseLoci)
                        varBeta[mSet][r] = sampleVarBeta(scaleM[pos],dfM[pos],beta[pos,theseLoci],regionSize)
                end
        end
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
