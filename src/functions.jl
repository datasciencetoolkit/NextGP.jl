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
export sampleBayesPR!,sampleBayesB!,sampleBayesC!,sampleBayesR!,sampleBayesLV!
export sampleZandZVar!

#Sampling fixed effects
function sampleX!(xSet::Union{Symbol,Tuple},X::Dict,b::Vector,ycorr::Vector,varE::Float64)
	if length(b[X[xSet].pos])==1
		ycorr    .+= X[xSet].data .* b[X[xSet].pos]
		rhs      = X[xSet].data'*ycorr .+ X[xSet].rhs
		meanMu   = X[xSet].ixpx*rhs			
                b[X[xSet].pos] .= rand(Normal(meanMu[],sqrt((X[xSet].ixpx*varE))[]))
		ycorr    .-= X[xSet].data .* b[X[xSet].pos]
	else
		ycorr    .+= X[xSet].data*b[X[xSet].pos]
                rhs      = X[xSet].data'*ycorr .+ X[xSet].rhs
                meanMu   = X[xSet].ixpx*rhs
		b[X[xSet].pos] .= rand(MvNormal(vec(meanMu),convert(Array,Symmetric(X[xSet].ixpx*varE))))
		ycorr    .-= X[xSet].data*b[X[xSet].pos]
	end
end

#sample random effects
function sampleU(zSet::Union{Expr,Symbol},Z::Dict,varE::Float64,varU::Dict,u::Vector,ycorr::Vector{Float64})
	uVec = deepcopy(u[Z[zSet].pos])
	λz = varE/varU[zSet]
	Yi = Z[zSet].Zp*ycorr #computation of Z'ycorr for ALL  rhsU
	nCol = length(uVec)
	for i in 1:nCol
        	uVec[i] = 0.0 #also excludes individual from iMat! Nice trick.
		rhsU = Yi[i] - λz*dot(view(Z[zSet].iVarStr,:,i),uVec)
                lhsU = getindex(Z[zSet].zpz,i) + (view(Z[zSet].iVarStr,i,i)*λz)[1]
		invLhsU = 1.0/lhsU
                meanU = invLhsU*rhsU
                uVec[i] = rand(Normal(meanU,sqrt(invLhsU*varE)))
        end
	return uVec
end


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


function sampleZ!(zSet::Union{Expr,Symbol},Z::Dict,u::Vector,ycorr::Vector{Float64},varE::Float64,varU::Dict)
        #for each random effect
	ycorr .+= Z[zSet].data*u[Z[zSet].pos]'
        u[Z[zSet].pos] .= sampleU(zSet,Z,varE,varU,u,ycorr)
	ycorr .-= Z[zSet].data*u[Z[zSet].pos]'		
	varU[zSet] = sampleVarU(Z[zSet].iVarStr,Z[zSet].scale,Z[zSet].df,u[Z[zSet].pos])		
end

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
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize::Int64 = length(theseLoci)
		lambda = varE/(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = BLAS.dot(view(M[mSet].data,:,locus),ycorr) + getindex(M[mSet].rhs,locus)
			lhs = getindex(M[mSet].mpm,locus) + lambda
			meanBeta = lhs\rhs
			setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
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
			RHS = ((getindex(M[mSet].Mp,locus)*ycorr)./varE) .+ getindex(M[mSet].rhs,locus)
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
		lambda = varE/(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = BLAS.dot(view(M[mSet].data,:,locus),ycorr) #+ getindex(M[mSet].rhs,locus)
			v0 = getindex(M[mSet].mpm,locus)*varE
			v1 = (getindex(M[mSet].mpm,locus)^2)*varBeta[mSet][r] + v0
        		logDelta0 = -0.5*(log(v0) + (rhs^2)/v0) + M[mSet].logPiOut            # this locus not fitted
			logDelta1 = -0.5*(log(v1) + (rhs^2)/v1) + M[mSet].logPiIn             # this locus fitted       
        		probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1))
			if rand() < probDelta1
				setindex!(delta[M[mSet].pos],1,locus)
				nLoci += 1
				lhs = getindex(M[mSet].mpm,locus) + lambda
				meanBeta = lhs\rhs
				setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
				@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),1)
			else 
				setindex!(beta[M[mSet].pos],0.0,locus)
				setindex!(delta[M[mSet].pos],0,locus)
				@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),0)
			end
		end
	end
#	println("pi=$(nLoci/M[mSet].dims[2])")
end

function sampleBayesC!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	nLoci = 0
	lambda = varE/(varBeta[mSet][1])
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesC, so r=locus
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = BLAS.dot(view(M[mSet].data,:,locus),ycorr) #+ getindex(M[mSet].rhs,locus)
			v0 = getindex(M[mSet].mpm,locus)*varE
			v1 = (getindex(M[mSet].mpm,locus)^2)*varBeta[mSet][1] + v0
        		logDelta0 = -0.5*(log(v0) + (rhs^2)/v0) + M[mSet].logPiOut            # this locus not fitted
			logDelta1 = -0.5*(log(v1) + (rhs^2)/v1) + M[mSet].logPiIn             # this locus fitted       
        		probDelta1 = 1.0/(1.0 + exp(logDelta0-logDelta1))
			if rand() < probDelta1
				setindex!(delta[M[mSet].pos],1,locus)
				nLoci += 1
				lhs = getindex(M[mSet].mpm,locus) + lambda
				meanBeta = lhs\rhs
				setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			else 
				setindex!(beta[M[mSet].pos],0.0,locus)
				setindex!(delta[M[mSet].pos],0,locus)
			end
		end
	end
	@inbounds varBeta[mSet][1] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,beta[M[mSet].pos],nLoci)
#	println("pi=$(nLoci/M[mSet].dims[2])")
end

function sampleBayesR!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local meanBeta::Float64
	nVarComp = length(M[mSet].vClass)
	nLoci = zeros(nVarComp)
	varc = varBeta[mSet][1].*M[mSet].vClass
	sumS = 0
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = BLAS.dot(view(M[mSet].data,:,locus),ycorr) #+ getindex(M[mSet].rhs,locus)
			lhs = zeros(nVarComp)
			ExpLogL = zeros(nVarComp)
			for v in 1:nVarComp
				lhs[v] = getindex(M[mSet].mpm,locus) + varE/varc[v]
				logLc = varc[v]==0.0 ? M[mSet].logPi[v] : -0.5*(log(varc[v]*lhs[v]/varE)-((rhs^2)/(varE*lhs[v]))) + M[mSet].logPi[v]
#				logLc = -0.5*(log(varc[v]*lhs[v]/varE)-((rhs^2)/(varE*lhs[v]))) + M[mSet].logPi[v]
				ExpLogL[v] = exp(logLc)
			end
			
			probs = ExpLogL./sum(ExpLogL)
			cumProbs = cumsum(probs)
			classSNP = findfirst(x->x>=rand(), cumProbs) #position
			setindex!(delta[M[mSet].pos],classSNP,locus)
			nLoci[classSNP] += 1
			meanBeta = lhs[classSNP]\rhs
			betaSample = sampleBeta(meanBeta, lhs[classSNP], varE)
			setindex!(beta[M[mSet].pos],betaSample,locus)
			BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
#			sumS += (betaSample^2)/M[mSet].vClass[classSNP]
		end
	end
	varSNP = getindex.(Ref(M[mSet].vClass),delta[M[mSet].pos][1,:])
	nonZeroPos = findall(!iszero, varSNP)
	nonZeroBeta = getindex.(Ref(beta[M[mSet].pos]),nonZeroPos)
	sumS = sum((nonZeroBeta.^2)./varSNP[nonZeroPos])
	
	@inbounds varBeta[mSet][1] = sampleVarBetaR(M[mSet].scale,M[mSet].df,sumS,length(nonZeroPos))
#	println("pi=$(nLoci./M[mSet].dims[2])")
#	println("var=$(varBeta[mSet][1].*M[mSet].vClass)")
end

function sampleBayesLV!(mSet::Symbol,M::Dict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	nLoci = 0
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		lambda = varE/(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = BLAS.dot(view(M[mSet].data,:,locus),ycorr)
			lhs = getindex(M[mSet].mpm,locus) + lambda
			meanBeta = lhs\rhs
			setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
			BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
		end
	end
	
	# model variance
	var_var = 0.001
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		for locus in theseLoci::UnitRange{Int64}
			vari = varBeta[mSet][locus]
			bi = getindex(beta[M[mSet].pos],locus)
			log_vari = log(vari)
			var_resid = M[mSet].SNPVARRESID[locus]
			var_mui = log_vari - var_resid
			
			c1 = ^(vari,-1.51)
			c2 = exp(-0.5*bi*bi/vari)*rand()
			c3 = exp(-0.5*var_resid*var_resid/var_var)*rand()
			temp = sqrt(-2*var_var*log(c3))
			lbound = exp(mu_var-temp)
			ubound = exp(mu_var+temp)
			(temp=exp((-2/3)*log(c1))) < rbound ? rbound=temp : nothing
			(temp= -0.5*bi*bi/log(c2)) > lbound ? lbound=temp : nothing
			if lboundd >= rbound
				println("Trap in sampling!!")
			else
				varBeta[mSet][locus] = lbound+rand()*(rbound-lbound)
				log_vari = log(varBeta[mSet][locus])
				M[mSet].SNPVARRESID[locus] = log_vari - mu_var
			end
			
			M[mSet].SNPVARRESID[locus] .+= M[mSet].covariates*c			
			rhsC = transpose(M[mSet].covariates)*ycorr
			meanC   = M[mSet].iCpC*rhsC
			c .= rand(MvNormal(vec(meanC),convert(Array,Symmetric(M[mSet].iCpC*var_var))))
			M[mSet].SNPVARRESID[locus] .-= M[mSet].covariates*c
		end
	end
end

#####



#sample marker effects
function sampleBeta(meanBeta, lhs, varE)
    return rand(Normal(meanBeta,sqrt(lhs\varE)))
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



end
