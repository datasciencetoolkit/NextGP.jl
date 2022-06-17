module samplers

using Distributions, LinearAlgebra
using StatsBase
using Printf
using CSV
using DataFrames
using DataStructures

include("outFiles.jl")
include("misc.jl")

export runSampler

#main sampler
function runSampler(rowID,A,Y,X,Z,chainLength,burnIn,outputFreq,priorVCV,M,paths2maps,rS)
	
	println("priorVCV $(priorVCV)")
	#output settings
	these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength) #print these iterations        

        #some info
	##This is not really nFix, but the "blocks" of fixed effects
        nFix  = length(X)
	nRand = length(Z)
	nData = length(Y)
	nMarkerSets = length(M)
	nMarkers    = [size(M[m],2) for m in keys(M)]
	println("nMarkers: $nMarkers")

        #initial computations and settings
	ycorr = deepcopy(Y)
	        
	##make iXpX, Z', zpz (for uncor)
        iXpX = deepcopy(X)
        for x in keys(X)
                iXpX[x] = inv(X[x]'X[x])
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
	println("XKeyPos: $XKeyPos")

	
	ZKeyPos = OrderedDict{Any,Int64}()
	[ZKeyPos[collect(keys(Z))[i]]=i for i in 1:length(keys(Z))]
        println("ZKeyPos: $ZKeyPos")


	
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
				iVarStr[zSet] = inv(A)
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
   	νS_E = scaleE*dfE



	#ADD MARKERS
	# read map file and make regions
	regionArray = OrderedDict{Any,Array{UnitRange{Int64},1}}()
        for mSet in keys(M)
                theseRegions = prep2RegionData(paths2maps[mSet],rS[mSet]) ###first data
                regionArray[mSet] = theseRegions
        end


        nRegions  = [length(regionArray[mSet]) for mSet in keys(regionArray)] #per component


	#make mpm
       	mpm = OrderedDict{Any,Any}()
#        for mSet in keys(M)
#		tempmpm = []
#		nowM = M[mSet]
#			for c in eachcol(nowM)
#				push!(tempmpm,BLAS.dot(c,c))
#			end
#		mpm[mSet] = tempmpm
#        end

#	println("mpm1: $(mpm["M1"][1:5]) $(mpm["M2"][1:5]) $(mpm["M3"][1:5])")

	##########
	
	corM = []

	for pSet ∈ keys(priorVCV)
		corEffects = []
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
			end
		else println("$pSet will be correlated")
			correlate = collect(pSet)
			for pSet in correlate
				println(pSet)
				push!(corEffects,pSet)
			end
			if issubset(corEffects,collect(keys(M)))
				push!(corM,pSet)
				mpm[pSet] = MatByMat.(hcat.(eachcol.(getindex.(Ref(M), (pSet)))...))
			end
		end
	end  	

	println("corM: $corM")

	dfM = Dict{String,Any}()	
	for mSet ∈ keys(mpm)
		dfM[mSet] = 3.0+size(priorVCV[mSet],1)
	end

	println("dfM $dfM")

	scaleM = Dict{Any,Any}()
        for mSet in keys(mpm)
                nMComp = size(priorVCV[mSet],1)
                nMComp > 1 ? scaleM[mSet] = priorVCV[mSet].*(dfM[mSet]-nMComp-1.0)  : scaleM[mSet] = priorVCV[mSet]*(dfM[mSet]-2.0)/dfM[mSet] #I make float and array of float
        end
	
	println("scaleM $scaleM")
	
	#########	

	#key positions for each effect, for speed. Order of matrices in M are preserved here.
	MKeyPos = OrderedDict{String,Int64}()
	for mSet in keys(M)
		pos = findall(mSet.==collect(keys(M)))[]
		MKeyPos[mSet] = pos
	end
	println("MKeyPos: $MKeyPos")	
	
	#storage

	varU = varU_prior #for storage


	beta = zeros(Float64,nMarkerSets,maximum(vcat([0,nMarkers]...))) #zero is for max to work when no SNP data is present #can allow unequal length! Remove tail zeros for printing....
#	vcovBeta = fill(Matrix(Diagonal(varM)),maximum(nRegions)) #can allow unequal length! Remove tail zeros for printing....


	varBeta = OrderedDict{Any,Any}()
	for mSet in keys(M)
		varBeta[mSet] = hcat(fill(priorVCV[mSet],nRegions[MKeyPos[mSet]])...) #later, direct reference to key when varM_prior is a dictionary
	end
	println("keys of varBeta: $(keys(varBeta))")

	#Start McMC
        for iter in 1:chainLength
		#sample residual variance
	       	varE = sampleVarE(νS_E,ycorr,dfE,nData)
		
		#sample fixed effects
        	#always returns corrected Y and new b
	        sampleX!(X,b,iXpX,nFix,nColEachX,XKeyPos,ycorr,varE)
	
		#sample random effects
		# always returns corrected Y and new u
		sampleZ!(iVarStr,Z,Zp,zpz,nRand,ZKeyPos,varE,varU,u,ycorr)

		#sample variances
		sampleRanVar!(varU,nRand,scaleU,dfDefault,u,ZKeyPos,iVarStr)
		
		#sample marker effects and variances
	        sampleMandMVar_view!(M,beta,mpm,nMarkerSets,MKeyPos,regionArray,nRegions,ycorr,varE,varBeta,scaleM,dfM)
               		
        	#print
		if iter in these2Keep
			IO.outMCMC(pwd(),"b",vcat(b...)') ### currently no path is provided!!!!
			IO.outMCMC(pwd(),"u",vcat(u...)')
			IO.outMCMC(pwd(),"varE",varE)
			IO.outMCMC(pwd(),"varU",hcat([varU[k] for k in keys(varU)]...))
			for markers in keys(M)
				IO.outMCMC(pwd(),"var"*markers,varBeta[markers])
			end
	#		if onScreen==true
            			println("b, $(vcat(b...))") #i always vectorize b. maybe better to make it vector initially

				for mSet in keys(M)
        	                	println("vG, $(var(M[mSet]*beta[MKeyPos[mSet],:]))")
				end


        #		end
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
		λz = varE/(varU[zSet])
	        ycorr .+= Zmat[zSet]*uVec		
	        Yi = ZpMat[zSet]*ycorr #computation of Z'ycorr for ALL  rhsU
		nCol = length(zpzMat[zSet])
	        for i in 1:nCol
        	        uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
			rhsU = Yi[i] - λz*dot(iMat[:,i],uVec)
                	lhsU = tempzpz[i] + (iMat[i,i]*λz)[1]
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

function sampleMandMVar_view!(MMat,beta,mpmMat,nMSet,keyM,regionsMat,regions,ycorr,varE,varBeta,scaleM,dfM)
        #for each marker set
        for mSet in keys(MMat)
		nowM = MMat[mSet]
                pos = keyM[mSet]
                for r in 1:regions[pos]
                        theseLoci = regionsMat[mSet][r]
                        regionSize = length(theseLoci)
                        lambda = varE/(varBeta[mSet][r])
                        for locus in theseLoci
                                BLAS.axpy!(beta[pos,locus],view(nowM,:,locus),ycorr)
                                rhs::Float64 = BLAS.dot(view(nowM,:,locus),ycorr)
                                lhs::Float64 = mpmMat[mSet][locus] + lambda
                                meanBeta::Float64 = lhs\rhs
                                beta[pos,locus] = sampleBeta(meanBeta, lhs, varE)
                                BLAS.axpy!(-1.0*beta[pos,locus],view(nowM,:,locus),ycorr)
                        end
                        varBeta[mSet][r] = sampleVarBeta(scaleM[mSet],dfM[mSet],beta[pos,theseLoci],regionSize)
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


#Sample residual variance
function sampleVarE(νS_e,yCorVec,df_e,nRecords)
    return((νS_e + BLAS.dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords)))
end



end
