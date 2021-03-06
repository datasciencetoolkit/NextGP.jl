module samplers


using Distributions, LinearAlgebra
using StatsBase
using Printf
using CSV
using DataFrames
using DataStructures
using ProgressMeter

include("outFiles.jl")
include("misc.jl")

export runSampler

#main sampler
function runSampler(rowID,iA,Y,X,Z,chainLength,burnIn,outputFreq,priorVCV,M,paths2maps,rS,outPut)
	
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	nRand = length(Z)
	nData = length(Y)
	nMarkerSets = length(M)
	nMarkers    = [size(M[m],2) for m in keys(M)]
	println("number of markers: $nMarkers")

        #initial computations and settings
	ycorr = deepcopy(Y)
	        
	##make iXpX, Z', zpz (for uncor)
        iXpX = deepcopy(X)
        for x in keys(X)
		XpX = X[x]'X[x]
#		if !isposdef(XpX)
			XpX += Matrix(I*minimum(abs.(diag(XpX)./size(X[x],1))),size(XpX))
#		end
               	iXpX[x] = inv(XpX)
        end

	Zp  = deepcopy(Z) #similar(Z')
#	zpz = Array{Array{Float64, 1},1}(undef,0)
	zpz = OrderedDict{Any,Any}()
	for z in keys(Z)
                zpz[z] = diag(Z[z]'Z[z])
		Zp[z]  = Z[z]'
        end
	
		
        #key positions for speed
        XKeyPos = OrderedDict{Any,Int64}()
        [XKeyPos[collect(keys(X))[i]]=i for i in 1:length(keys(X))]

	
	ZKeyPos = OrderedDict{Any,Int64}()
	[ZKeyPos[collect(keys(Z))[i]]=i for i in 1:length(keys(Z))]
	
        ##make b and u arrays
        b = Array{Array{Float64, 1},1}(undef,0)
        ##counts columns per effect
        nColEachX = []
        for xSet in keys(X)
                println(xSet)
                nCol = size(X[xSet],2)
                push!(b,fill(0.0,nCol))
                nColEachX = push!(nColEachX,nCol)
        end

        u = Array{Array{Float64, 1},1}(undef,0)
        ##counts columns per effect
        nColEachZ = []
	##get priors per effect
	iVarStr = Dict{Any,Array{Float64,2}}() #inverses will be computed
	varU_prior = OrderedDict{Any,Any}()
        for zSet in keys(Z)
                nCol = size(Z[zSet],2)
                push!(u,fill(0.0,nCol))
                nColEachZ = push!(nColEachZ,nCol)
		#var structures and priors
		if haskey(priorVCV,zSet)	
			if isempty(priorVCV[zSet][1])
				println("priorVCV structure for $zSet is empty, an identity matrix will be used")
				iVarStr[zSet] = Matrix(1.0I,nCol,nCol)
			elseif priorVCV[zSet][1]=="A"
				iVarStr[zSet] = iA
				println("priorVCV structure for $zSet is A, computed A matrix will be used")
			else 	iVarStr[zSet] = inv(priorVCV[zSet][1])
			end
			varU_prior[zSet] = priorVCV[zSet][2]
		else	println("priorVCV structure for $zSet is empty, an identity matrix will be used with an arbitrary variance of 100")
			iVarStr[zSet] = Matrix(1.0I,nCol,nCol)
			varU_prior[zSet] = 100	
		end
        end

	#set up for E	
	isempty(priorVCV["e"][1]) ? strE = Matrix(1.0I,nData,nData) : strE = priorVCV["e"][1]
	varE_prior = priorVCV["e"][2]

	#parameters for priors
        dfE = 4.0
	dfDefault = 4.0
 
	       
	if varE_prior==0.0
		varE_prior  = 0.0005
       		scaleE     = 0.0005
        else
       		scaleE    = varE_prior*(dfE-2.0)/dfE    
   	end
	
	##no correlated random effects
	scaleU = Dict{Any,Any}()
	for zSet in keys(Z)
		scaleU[zSet] = varU_prior[zSet]*(dfDefault-2.0)/dfDefault
	end	


	#pre-computations using priors, not relevant for correlated random effects
   	??S_E = scaleE*dfE



	#ADD MARKERS
	# read map file and make regions

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

	for pSet ??? keys(priorVCV)
		corEffects = []
		corPositions = []
		if typeof(pSet)==String
			println("$pSet is univariate")
			if pSet ??? keys(M)
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

#	nRegions  = [mSet => length(regionArray[mSet]) for mSet in keys(regionArray)] #per marker set
	nRegions  = OrderedDict(mSet => length(regionArray[mSet]) for mSet in keys(regionArray))
	println("number of regions: $(values(nRegions))")


	dfM = Dict{Any,Any}()	
	for mSet ??? keys(mpm)
		dfM[mSet] = 3.0+size(priorVCV[mSet],1)
	end


	scaleM = Dict{Any,Any}()
        for mSet in keys(mpm)
                nMComp = size(priorVCV[mSet],1)
                nMComp > 1 ? scaleM[mSet] = priorVCV[mSet].*(dfM[mSet]-nMComp-1.0)  : scaleM[mSet] = priorVCV[mSet]*(dfM[mSet]-2.0)/dfM[mSet] #I make float and array of float
        end
	
	
	#storage

	varU = varU_prior #for storage


	beta = zeros(Float64,nMarkerSets,maximum(vcat([0,nMarkers]...))) #zero is for max to work when no SNP data is present #can allow unequal length! Remove tail zeros for printing....
#	vcovBeta = fill(Matrix(Diagonal(varM)),maximum(nRegions)) #can allow unequal length! Remove tail zeros for printing....


        varBeta = OrderedDict{Any,Any}()
        for mSet in keys(mpm)
                varBeta[mSet] = [priorVCV[mSet] for i in 1:length(regionArray[mSet])] #later, direct reference to key when varM_prior is a dictionary
        end


	#Start McMC
@showprogress 1 "MCMC progress..." for iter in 1:chainLength
	sleep(0.1)
	
		#sample residual variance
	       	varE = sampleVarE(??S_E,ycorr,dfE,nData)
		
		#sample fixed effects
        	#always returns corrected Y and new b
	        sampleX!(X,b,iXpX,nFix,nColEachX,XKeyPos,ycorr,varE)
	
		#sample random effects
		# always returns corrected Y and new u
		sampleZ!(iVarStr,Z,Zp,zpz,nRand,ZKeyPos,varE,varU,u,ycorr)

		#sample variances
		sampleRanVar!(varU,nRand,scaleU,dfDefault,u,ZKeyPos,iVarStr)
		
		#sample marker effects and variances
	        sampleMandMVar_view!(M,Mp,corM,corMPos,beta,mpm,nMarkerSets,BetaKeyPos,regionArray,nRegions,ycorr,varE,varBeta,scaleM,dfM)
               		
        	#print
		if iter in these2Keep
			IO.outMCMC(outPut,"b",vcat(b...)') ### currently no path is provided!!!!
			IO.outMCMC(outPut,"u",vcat(u...)')
			IO.outMCMC(outPut,"varE",varE)
			IO.outMCMC(outPut,"varU",hcat([varU[k] for k in keys(varU)]...))
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
                ycorr    .+= X[xSet]*b[pos]
                rhs      = X[xSet]'*ycorr
                meanMu   = iXpX[xSet]*rhs
                if nColEachX[pos] == 1
                        b[pos] .= rand(Normal(meanMu[],sqrt((iXpX[xSet]*varE))[]))
                else b[pos] .= rand(MvNormal(vec(meanMu),convert(Array,Symmetric(iXpX[xSet]*varE))))
                end
                ycorr    .-= X[xSet]*b[pos]
        end
end

#Sampling random effects
function sampleZ!(iStrMat,Zmat,ZpMat,zpzMat,nRand,ZKeyPos,varE,varU,u,ycorr)
	#block for each effect
	for zSet in keys(Zmat)
		pos = ZKeyPos[zSet]
		uVec = deepcopy(u[pos])
		iMat = iStrMat[zSet]
		tempzpz = zpzMat[zSet] ###added
		??z = varE/(varU[zSet])
	        ycorr .+= Zmat[zSet]*uVec		
	        Yi = ZpMat[zSet]*ycorr #computation of Z'ycorr for ALL  rhsU
		nCol = length(zpzMat[zSet])
	        for i in 1:nCol
        	        uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
			rhsU = Yi[i] - ??z*dot(iMat[:,i],uVec)
                	lhsU = tempzpz[i] + (iMat[i,i]*??z)[1]
			invLhsU = 1.0/lhsU
                	meanU = invLhsU*rhsU
                	uVec[i] = rand(Normal(meanU,sqrt(invLhsU*varE)))
        	end
		u[pos] = uVec
        	ycorr .-= Zmat[zSet]*uVec
	end
end


#Sampling marker effects
function sampleM!(MMat,beta,mpmMat,nMSet,keyM,regionsMat,regions,ycorr,varE,varBeta)
        #for each marker set
        for mSet in keys(MMat)
                for r in 1:regions[keyM[mSet]]
                        theseLoci = regionsMat[keyM[mSet]][r]
                        regionSize = length(theseLoci)
                        lambda = varE/(varBeta[keyM[mSet]][r])
                        for locus in theseLoci
                                BLAS.axpy!(beta[keyM[mSet],locus],MMat[mSet][:,locus],ycorr)
                                rhs = BLAS.dot(MMat[mSet][:,locus],ycorr)
                                lhs = mpmMat[mSet][locus] + lambda
                                meanBeta = lhs\rhs
                                beta[keyM[mSet],locus] = sampleBeta(meanBeta, lhs, varE)
                                BLAS.axpy!(-1.0*beta[keyM[mSet],locus],MMat[mSet][:,locus],ycorr)
                        end
                end
        end
end

function sampleMandMVar!(MMat,beta,mpmMat,nMSet,keyM,regionsMat,regions,ycorr,varE,varBeta,scaleM,dfM)
        #for each marker set
        for mSet in keys(MMat)
		pos = keyM[mSet]
                for r in 1:regions[pos]
                        theseLoci = regionsMat[pos][r]
                        regionSize = length(theseLoci)
                        lambda = varE/(varBeta[mSet][r])
                        for locus in theseLoci
                                BLAS.axpy!(beta[pos,locus],MMat[mSet][:,locus],ycorr)
                                rhs = BLAS.dot(MMat[mSet][:,locus],ycorr)
                                lhs = mpmMat[mSet][locus] + lambda
                                meanBeta = lhs\rhs
                                beta[pos,locus] = sampleBeta(meanBeta, lhs, varE)
                                BLAS.axpy!(-1.0*beta[pos,locus],MMat[mSet][:,locus],ycorr)
                        end
			varBeta[mSet][r] = sampleVarBeta(scaleM[pos],dfM[pos],beta[pos,theseLoci],regionSize)
                end
        end
end

function sampleMandMVar_view!(MMat,MpMat,correlatedM,keyCorM,beta,mpmMat,nMSet,keyBeta,regionsMat,regions,ycorr,varE,varBeta,scaleM,dfM)
        #for each marker set
        for mSet in keys(mpmMat)
		if mSet in keys(correlatedM)
#			println("$mSet in corM")
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

#sample random effects' variances
function sampleRanVar!(varU,nRand,scale_ranVar,df_ranVar,effVec,keyZ,iStrMat)
	for z in keys(scale_ranVar)
		pos = keyZ[z]
		n = size(iStrMat[z],2)
		varU[z] = (scale_ranVar[z]*df_ranVar + effVec[pos]'*iStrMat[z]*effVec[pos])/rand(Chisq(df_ranVar + n))
	end
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
function sampleVarE(??S_e,yCorVec,df_e,nRecords)
    return (??S_e + BLAS.dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords))
end



end
