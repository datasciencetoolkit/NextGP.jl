module functions


using Distributions, LinearAlgebra
using StatsBase
using Printf
using CSV
using DataFrames
using DataStructures

include("outFiles.jl")

export sampleVarE
export sampleX!
export sampleBayesPR!,sampleBayesB!,sampleBayesC!,sampleBayesR!,sampleBayesRCπ!,sampleBayesLV!,sampleBayesRCplus!
export sampleZandZVar!

#Sampling fixed effects

### NEW, Wang's trick

function sampleb!(xSet::Union{Symbol,Tuple},X::Dict,b::Vector,ycorr::Vector,varE::Float64)
	iVarE = inv(varE)
	bVec = deepcopy(b[X[xSet].pos])
	Yi = X[xSet].Xp*ycorr*iVarE #computation of X'ycorr*iVarE for ALL  rhsb
	nCol = length(bVec)
	for i in 1:nCol
        	bVec[i] = 0.0 #also excludes the effect from iMat! Nice trick.
		rhsb = Yi[i] - dot(view(X[xSet].xpx,i,:),bVec)*iVarE
                lhsb = getindex(X[xSet].xpx,i,i)*iVarE
		invLhsb = 1.0/lhsb
                meanb = invLhsb*rhsb
                bVec[i] = rand(Normal(meanb,sqrt(invLhsb)))
        end
	return bVec
end

# NEW with D and with Wang's Trick
function sampleX!(xSet::Union{Symbol,Tuple},X::Dict,b::Vector,ycorr::Vector,varE::Float64)
	iVarE = inv(varE)
	if length(b[X[xSet].pos])==1
		ycorr    .+= X[xSet].data .* b[X[xSet].pos]
		rhs      = (X[xSet].Xp*ycorr).*iVarE .+ X[xSet].rhs
		lhs      = X[xSet].xpx .*iVarE .+ X[xSet].lhs
		meanMu   = lhs\rhs			
                b[X[xSet].pos] .= rand(Normal(meanMu[],sqrt(inv(lhs[]))))
		ycorr    .-= X[xSet].data .* b[X[xSet].pos]
	else
		ycorr    .+= X[xSet].data*b[X[xSet].pos]
		b[X[xSet].pos] .= sampleb!(xSet,X,b,ycorr,varE)
		ycorr    .-= X[xSet].data*b[X[xSet].pos]
	end
end

#sample random effects
#Uni D u
function sampleU(zSet::Union{Expr,Symbol},Z::Dict,varE::Float64,varU::Dict,u::Vector,ycorr::Vector{Float64})
	uVec = deepcopy(u[Z[zSet].pos])
	iVarE = 1/varE
	iVarU = 1/varU[zSet]
	Yi = Z[zSet].Zp*ycorr*iVarE #computation of Z'*D^-1*ycorr*iVarE for ALL  rhsU
	nCol = length(uVec)
	for i in 1:nCol
        	uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
		rhsU = Yi[i] - iVarU*dot(view(Z[zSet].iVarStr,:,i),uVec)
                lhsU = getindex(Z[zSet].zpz,i)*iVarE + (view(Z[zSet].iVarStr,i,i)*iVarU)[1]
		invLhsU = 1.0/lhsU
                meanU = invLhsU*rhsU
                uVec[i] = rand(Normal(meanU,sqrt(invLhsU)))
        end
	return uVec
end

#Mul u no D
function sampleU(zSet::Tuple,Z::Dict,varE::Float64,varU::Dict,u::Vector,ycorr::Vector{Float64})
	uVec = deepcopy(u[Z[zSet].pos])
	nCol = size(uVec,2)
	iVarU = inv(varU[zSet])
	for i in 1:nCol
		setindex!(uVec,[0;0],:,i)
		Yi = Z[zSet].Zp[i]*ycorr		
		rhsU = (Yi./varE) - kron(view(Z[zSet].iVarStr,[i],:),iVarU)*vec(uVec)
                invLhsU = inv((getindex(Z[zSet].zpz,i)./varE) + (view(Z[zSet].iVarStr,i,i).*iVarU))
                meanU = invLhsU*rhsU
		setindex!(uVec,rand(MvNormal(meanU,convert(Array,Symmetric(invLhsU)))),:,i)
        end
	return uVec
end

#Uni Main
function sampleZ!(zSet::Union{Expr,Symbol},Z::Dict,u::Vector,ycorr::Vector{Float64},varE::Float64,varU::Dict)
        #for each random effect
	ycorr .+= Z[zSet].data*u[Z[zSet].pos]'
        u[Z[zSet].pos] .= sampleU(zSet,Z,varE,varU,u,ycorr)
	ycorr .-= Z[zSet].data*u[Z[zSet].pos]'		
	varU[zSet] = sampleVarU(Z[zSet].iVarStr,Z[zSet].scale,Z[zSet].df,u[Z[zSet].pos])		
end

#Mul Main no D
function sampleZ!(zSet::Tuple,Z::Dict,u::Vector,ycorr::Vector{Float64},varE::Float64,varU::Dict)
	nCol = size(u[Z[zSet].pos],2)
	for i in 1:nCol
		ycorr .+= Z[zSet].data[i]*getindex(u[Z[zSet].pos],:,i)
	end
	u[Z[zSet].pos] .= sampleU(zSet,Z,varE,varU,u,ycorr)
	varU[zSet] = sampleCoVarU(Z[zSet].iVarStr,Z[zSet].scale,Z[zSet].df,u[Z[zSet].pos])
	for i in 1:nCol
		ycorr .-= Z[zSet].data[i]*getindex(u[Z[zSet].pos],:,i)
	end
end



#sample random marker effects

##### Component-wise, seperated functions for symbol and tuple

function sampleBayesPR!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	iVarE = 1/varE
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize::Int64 = length(theseLoci)
		iVarBeta = 1/varBeta[mSet][r]
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = getindex(M[mSet].Mp,locus)*ycorr.*iVarE + getindex(M[mSet].rhs,locus)
			lhs = getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + iVarBeta
			meanBeta = lhs\rhs
			setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs),locus)
			BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
		end
		@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),regionSize)
	end
end
	
##### Component-wise, seperated functions for symbol and tuple
function sampleBayesPR!(mSet::Tuple,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize = length(theseLoci)
		invB = inv(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			ycorr .+= M[mSet].data[locus]*getindex.(beta[M[mSet].pos],locus)
			RHS = ((getindex(M[mSet].Mp,locus)*ycorr)./varE)
			invLHS::Array{Float64,2} = inv((getindex(M[mSet].mpm,locus)./varE) .+ invB)
			meanBETA::Array{Float64,1} = invLHS*RHS
			setindex!.(beta[M[mSet].pos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),locus)
			ycorr .-= M[mSet].data[locus]*getindex.(beta[M[mSet].pos],locus)	
		end
		@inbounds varBeta[mSet][r] = sampleVarCovBetaPR(M[mSet].scale,M[mSet].df,reduce(hcat,getindex.(beta[M[mSet].pos],Ref(theseLoci))),regionSize)
	end
end


function sampleBayesB!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	nLoci = 0
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		iVarE = 1/varE
		iVarBeta = 1/varBeta[mSet][r]
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rrr = BLAS.dot(view(M[mSet].data,:,locus),ycorr) #not rhs!
			v0 = getindex(M[mSet].mpm,locus)*varE
			v1 = (getindex(M[mSet].mpm,locus)^2)*varBeta[mSet][r] + v0
        		logDelta0 = -0.5*(log(v0) + (rrr^2)/v0) + M[mSet].logPi[1]            # this locus not fitted
			logDelta1 = -0.5*(log(v1) + (rrr^2)/v1) + M[mSet].logPi[2]             # this locus fitted       
        		probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1))
			if rand() < probDelta1
				setindex!(delta[M[mSet].pos],1,locus)
				nLoci += 1
				rhs = getindex(M[mSet].Mp,locus)*ycorr*iVarE + getindex(M[mSet].rhs,locus)
				lhs = getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + iVarBeta
				meanBeta = lhs\rhs
				setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs),locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
				@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),1)
			else 
				setindex!(beta[M[mSet].pos],0.0,locus)
				setindex!(delta[M[mSet].pos],0,locus)
				@inbounds varBeta[mSet][r] = 0.0
			end
		end
	end
	if M[mSet].estPi == true 
		piIn = samplePi(nLoci,M[mSet].dims[2]) #probability of in
		M[mSet].piHat .= [1.0-piIn piIn]
		M[mSet].logPi .= log.([1.0-piIn piIn])
	end
end

function sampleBayesC!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	nLoci = 0
	iVarE = 1/varE
	iVarBeta = 1/varBeta[mSet][1]
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesC, so r=locus
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rrr = BLAS.dot(view(M[mSet].data,:,locus),ycorr)
			v0 = getindex(M[mSet].mpm,locus)*varE
			v1 = (getindex(M[mSet].mpm,locus)^2)*varBeta[mSet][1] + v0

			logDelta0 = -0.5*(log(v0) + (rrr^2)/v0) + M[mSet].logPi[1]            # this locus not fitted
			logDelta1 = -0.5*(log(v1) + (rrr^2)/v1) + M[mSet].logPi[2]             # this locus fitted       

        		probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1))
			if rand() < probDelta1
				setindex!(delta[M[mSet].pos],1,locus)
				nLoci += 1
				rhs = getindex(M[mSet].Mp,locus)*ycorr*iVarE #+ getindex(M[mSet].rhs,locus)
				lhs = getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + iVarBeta
				meanBeta = lhs\rhs
				setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs),locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			else 
				setindex!(beta[M[mSet].pos],0.0,locus)
				setindex!(delta[M[mSet].pos],0,locus)
			end
		end
	end
	@inbounds varBeta[mSet][1] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,beta[M[mSet].pos],nLoci)
	if M[mSet].estPi == true 
		piIn = samplePi(nLoci,M[mSet].dims[2]) #probability of in
		M[mSet].piHat .= [1.0-piIn piIn]
		M[mSet].logPi .= log.([1.0-piIn piIn])
	end
end

function sampleBayesR!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local meanBeta::Float64
	nVarClass = length(M[mSet].vClass)
	nLoci     = zeros(Int64,nVarClass)
	nNonZero = 0
	varc      = varBeta[mSet][1].*M[mSet].vClass
	sumS      = 0
	iVarE = 1/varE
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesR
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = getindex(M[mSet].Mp,locus)*ycorr*iVarE + getindex(M[mSet].rhs,locus)
			lhs = zeros(nVarClass)
			ExpLogL = zeros(nVarClass)
			for v in 1:nVarClass
				lhs[v] = varc[v]==0.0 ? 0.0 : getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + 1/varc[v]
				logLc  = varc[v]==0.0 ? M[mSet].logPi[v] : -0.5*(log(varc[v]*lhs[v])-((rhs^2)/lhs[v])) + M[mSet].logPi[v]
				ExpLogL[v] = exp(logLc)
			end
			
			probs = ExpLogL./sum(ExpLogL)
			cumProbs = cumsum(probs)
			classSNP = findfirst(x->x>=rand(), cumProbs) #position
			setindex!(delta[M[mSet].pos],classSNP,locus)
			nLoci[classSNP] += 1
			###sample only non-zero class SNPs
			if varc[classSNP]!= 0.0
				nNonZero += 1
				meanBeta = lhs[classSNP]\rhs
				betaSample = sampleBeta(meanBeta, lhs[classSNP])
				setindex!(beta[M[mSet].pos],betaSample,locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
				##
				varSNP = M[mSet].vClass[classSNP]
				sumS +=  betaSample^2 / varSNP  
				##
			else setindex!(beta[M[mSet].pos],0.0,locus)
			end
		end
	end
		
	##
	@inbounds varBeta[mSet][1] = sampleVarBetaR(M[mSet].scale,M[mSet].df,sumS,nNonZero)
	##
	
	if M[mSet].estPi == true 
		piHat = samplePi(nLoci)
		M[mSet].piHat .= piHat
		M[mSet].logPi .= log.(piHat)
	end
end

function sampleBayesRCπ!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local meanBeta::Float64
	nAnnot    = nVarCov = M[mSet].nVarCov
	nVarClass = length(M[mSet].vClass)
	nLoci     = zeros(Int64,nAnnot,nVarClass)
	nNonZero = zeros(Int64,nAnnot)
	varc      = [v.*M[mSet].vClass for v in varBeta[mSet]]
	sumS 	  = zeros(Float64,nAnnot)
	iVarE = 1/varE
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesR
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = getindex(M[mSet].Mp,locus)*ycorr*iVarE + getindex(M[mSet].rhs,locus)
			lhs = zeros(nAnnot,nVarClass)
			ExpLogL = zeros(nAnnot,nVarClass)
			for a in M[mSet].annotNonZeroPos[locus]
				for v in 1:nVarClass
					lhs[a,v] = varc[a][v]==0.0 ? 0.0 : getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + 1/varc[a][v]
					logLv    = varc[a][v]==0.0 ? M[mSet].logPi[a][v] : -0.5*(log(varc[a][v]*lhs[a,v])-((rhs^2)/lhs[a,v])) + M[mSet].logPi[a][v]
					ExpLogL[a,v] = exp(logLv)
				end
			end
			
			probAnnot1 = M[mSet].annotProb[locus,:] .* vec(sum(ExpLogL,dims=2))
			probAnnot2 = sum(probAnnot1)
			probAnnot = probAnnot1 ./ probAnnot2
			##########
			AnnnotClassSNP = rand(Categorical(probAnnot))  #position
			posAnnotInNonZero = findfirst(isequal(AnnnotClassSNP), M[mSet].annotNonZeroPos[locus])
			##pi sampled here
			M[mSet].annotProb[locus,M[mSet].annotNonZeroPos[locus]] = sampleProb(posAnnotInNonZero,M[mSet].annotInput[locus,M[mSet].annotNonZeroPos[locus]])
			##########
				
			probsV = ExpLogL[AnnnotClassSNP,:]./sum(ExpLogL[AnnnotClassSNP,:])
			cumProbsV = cumsum(probsV)
			classSNP = findfirst(x->x>=rand(), cumProbsV) #position
			
			setindex!(delta[M[mSet].pos],classSNP,locus)
			setindex!(M[mSet].annotCat,AnnnotClassSNP,locus)
			nLoci[AnnnotClassSNP,classSNP] += 1
			###sample only non-zero class SNPs
			if varc[AnnnotClassSNP][classSNP]!= 0.0
				nNonZero[AnnnotClassSNP] += 1
				meanBeta = lhs[AnnnotClassSNP,classSNP]\rhs
				betaSample = sampleBeta(meanBeta, lhs[AnnnotClassSNP,classSNP])
				setindex!(beta[M[mSet].pos],betaSample,locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
				varSNP = M[mSet].vClass[classSNP] #Same variant classes for all annotations
				sumS[AnnnotClassSNP] +=  betaSample^2 / varSNP  
			else setindex!(beta[M[mSet].pos],0.0,locus)
			end
		end
	end
		
	## Assumes same variant classes, v for all!
	for a in 1:nAnnot
		@inbounds varBeta[mSet][a] = sampleVarBetaR(M[mSet].scale,M[mSet].df,sumS[a],nNonZero[a])
	end
	
	#Estimate both Prob (p) and pi π
	if M[mSet].estPi == true
		for a in 1:nAnnot
			piHat = samplePi(nLoci[a,:])
			M[mSet].piHat[a] = piHat
			M[mSet].logPi[a] = log.(piHat)
		end
#		sampleProb()
	end
end

function sampleBayesRCplus!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local meanBeta::Float64
	nAnnot    = nVarCov = M[mSet].nVarCov
	nVarClass = length(M[mSet].vClass)
	nLoci     = zeros(Int64,nAnnot,nVarClass)
	nNonZero = zeros(Int64,nAnnot)
	varc      = [v.*M[mSet].vClass for v in varBeta[mSet]]
	sumS 	  = zeros(Float64,nAnnot)
	iVarE = 1/varE
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesR
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			lhs = zeros(nAnnot,nVarClass)
			ExpLogL = zeros(nAnnot,nVarClass)
			tempBeta = 0.0
			for a in M[mSet].annotNonZeroPos[locus]
				rhs = getindex(M[mSet].Mp,locus)*ycorr*iVarE + getindex(M[mSet].rhs,locus)
				for v in 1:nVarClass
					lhs[a,v] = varc[a][v]==0.0 ? 0.0 : getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + 1/varc[a][v]
					logLv    = varc[a][v]==0.0 ? M[mSet].logPi[a][v] : -0.5*(log(varc[a][v]*lhs[a,v])-((rhs^2)/lhs[a,v])) + M[mSet].logPi[a][v]
					ExpLogL[a,v] = exp(logLv)
				end
				probsV = ExpLogL[a,:]./sum(ExpLogL[a,:])
				cumProbsV = cumsum(probsV)
				classSNP = findfirst(x->x>=rand(), cumProbsV) #position
				setindex!(delta[M[mSet].pos],classSNP,locus)
				nLoci[a,classSNP] += 1
				###sample only non-zero class SNPs
				if varc[a][classSNP]!= 0.0
					nNonZero[a] += 1
					meanBeta = lhs[a,classSNP]\rhs
					betaSample = sampleBeta(meanBeta, lhs[a,classSNP])
					varSNP = M[mSet].vClass[classSNP] #Same variant classes for all annotations
					sumS[a] +=  betaSample^2 / varSNP  
				else betaSample = 0.0
				end
				tempBeta += betaSample
				BLAS.axpy!(-1.0*betaSample,view(M[mSet].data,:,locus),ycorr)
			end
			setindex!(beta[M[mSet].pos],tempBeta,locus)			
		end
	end
		
	## Assumes same variant classes, v for all!
	for a in 1:nAnnot
		@inbounds varBeta[mSet][a] = sampleVarBetaR(M[mSet].scale,M[mSet].df,sumS[a],nNonZero[a])
	end
	
	#Estimate both Prob (p) and pi π
	if M[mSet].estPi == true
		for a in 1:nAnnot
			piHat = samplePi(nLoci[a,:])
			M[mSet].piHat[a] = piHat
			M[mSet].logPi[a] = log.(piHat)
		end
	end
end

function sampleBayesLV!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64

	
#	var_var = sampleVarE(0,0.0,M[mSet].SNPVARRESID,length(M[mSet].SNPVARRESID)) #var(M[mSet].SNPVARRESID) #0.01
	var_var = M[mSet].estVarZeta == true ? var(M[mSet].SNPVARRESID) : M[mSet].varZeta[]
	println("var: $(var(beta[M[mSet].pos]))")

	setindex!(M[mSet].varZeta,var_var,1)
	
	nLoci = 0
	iVarE = 1/varE
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = getindex(M[mSet].Mp,locus)*ycorr*iVarE + getindex(M[mSet].rhs,locus)  
			lhs = getindex(M[mSet].mpm,locus)*iVarE + getindex(M[mSet].lhs,locus) + 1/varBeta[mSet][locus]
			meanBeta = lhs\rhs
			setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs),locus)
			BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
		end
	end
	
	# model variance
	
	trapped = 0
	notTrapped = 0
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		for locus in theseLoci::UnitRange{Int64}			
			vari = varBeta[mSet][locus]
			bi = getindex(beta[M[mSet].pos],locus)
			log_vari = M[mSet].logVar[locus]
			ζ = M[mSet].SNPVARRESID[locus]	#residual of variance for log-var
			var_mui = log_vari - ζ 		#mean of "variance at log scale"
			c1 = ^(vari,-1.5)*rand()
			c2 = exp(-0.5*bi*bi/vari)*rand()
			c3 = exp(-0.5*ζ*ζ/var_var)*rand()
			temp = sqrt(-2*var_var*log(c3))
			lbound = exp(var_mui-temp)
			rbound = exp(var_mui+temp)
			exp((-2/3)*log(c1)) < rbound ? rbound=exp((-2/3)*log(c1)) : nothing
			-0.5*bi*bi/log(c2) > lbound ? lbound=-0.5*bi*bi/log(c2) : nothing
			if lbound >= rbound
				trapped +=1
			else
				notTrapped +=1
				vari = lbound+rand()*(rbound-lbound)				
				varBeta[mSet][locus] = vari
				M[mSet].logVar[locus] = log(vari)
			end
		end
	end
	println("trapped: $(trapped/(trapped+notTrapped))")

	#rhsC = M[mSet].covariatesT*log.(varBeta[mSet])
	rhsC = M[mSet].covariatesT*M[mSet].logVar
	meanC = M[mSet].iCpC*rhsC
	M[mSet].c .= rand(MvNormal(vec(meanC),convert(Array,Symmetric(M[mSet].iCpC*var_var))))
	#M[mSet].SNPVARRESID .= log.(varBeta[mSet]) .- M[mSet].covariates*M[mSet].c
	M[mSet].SNPVARRESID .= M[mSet].logVar .- M[mSet].covariates*M[mSet].c
end

#####



#sample marker effects
function sampleBeta(meanBeta, lhs)
    return rand(Normal(meanBeta,sqrt(1/lhs)))
end

#sample random effects' variances (new U)
function sampleVarU(iMat,scale_ranVar,df_ranVar,effVec)
	n = size(iMat,2)
	return (scale_ranVar*df_ranVar + (effVec*iMat*effVec')[])/rand(Chisq(df_ranVar + n))
end

function sampleCoVarU(iMat,scale_ranVar,df_ranVar,effVec)
	n = size(iMat,2)
        return rand(InverseWishart(df_ranVar + n, convert(Array,Symmetric(effVec*iMat*effVec' + scale_ranVar))))
end

#sample marker variances
function sampleVarBetaPR(scalem,dfm,whichLoci,regionSize)::Float64
	return (scalem*dfm + BLAS.dot(whichLoci,whichLoci)) / rand(Chisq(dfm + regionSize))
end

function sampleVarCovBetaPR(scalem,dfm,whichLoci,regionSize)
	Sb = whichLoci'whichLoci
	return rand(InverseWishart(dfm + regionSize, scalem + Sb))
end
				
function sampleVarBetaR(scalem,dfm,sumS,nLoci)::Float64
	return (scalem*dfm + sumS) / rand(Chisq(dfm + nLoci))
end

#Sample residual variance
function sampleVarE(df_e,S_e,yCorVec,nRecords)
	return (df_e*S_e + BLAS.dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords))
end
function sampleVarE(E::NamedTuple,yCorVec,nRecords)
	return (E.df*E.scale + sum(E.iVarStr.*(yCorVec.^2)))/rand(Chisq(E.df + nRecords))
end
					
# +1 is for beta(1,1) prior
function samplePi(nIn::Int, nTotal::Int)
	return rand(Beta(nIn+1,nTotal-nIn+1))
end

# +1 is for Dirichlet(1,1,...,) prior						
function samplePi(nSNPs::Vector{Int})
	return rand(Dirichlet(nSNPs.+1))
end

# +1 already comes from the inputProb, posAnnotInNonZero comes from random sampling
function sampleProb(posAnnotInNonZero,input)
	input[posAnnotInNonZero]+=1
	return rand(Dirichlet(input))
end

							
end
