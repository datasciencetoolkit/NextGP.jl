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
export sampleBayesPR!
export BayesPR
export sampleZandZVar!

#Sampling fixed effects
function sampleX!(xMat,b,ycorr,varE)
	if length(b[xMat.pos])==1
		ycorr    .+= xMat.data .* b[xMat.pos]
		rhs      = xMat.data'*ycorr .+ xMat.rhs
		meanMu   = xMat.ixpx*rhs			
                b[xMat.pos] .= rand(Normal(meanMu[],sqrt((xMat.ixpx*varE))[]))
		ycorr    .-= xMat.data .* b[xMat.pos]
	else
		ycorr    .+= xMat.data*b[xMat.pos]
                rhs      = xMat.data'*ycorr .+ xMat.rhs
                meanMu   = xMat.ixpx*rhs
		b[xMat.pos] .= rand(MvNormal(vec(meanMu),convert(Array,Symmetric(xMat.ixpx*varE))))
		ycorr    .-= xMat.data*b[xMat.pos]
	end
end

#sample random effects
function sampleU(zSet::Union{Expr,Symbol},iMat,pos,ZComp,ZpComp,zpzComp,varE,varUComp,uVector,ycorr)
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


#sample random effects
function sampleU(zSet::Tuple,iMat,pos,ZComp,ZpComp,zpzComp,varE,varUComp,uVector,ycorr)
        uVec = deepcopy(uVector)
        λz = varE./varUComp
        Yi = ZpComp*ycorr #computation of Z'ycorr for ALL  rhsU
        nCol = length(zpzComp)
	for i in 1:nCol
                uVec[:,i] .= 0.0 #also excludes individual from iMat! Nice trick.
                rhsU = Yi[:,i] .- λz*dot(iMat[:,i],uVec)
                lhsU = zpzComp[i] + (iMat[i,i]*λz)[1]
                invLhsU = inv(lhsU)
                meanU = invLhsU*rhsU
                uVec[:,i] = rand(MvNormal(meanU,invLhsU*varE))
        end
        return uVec
end


function sampleZandZVar!(iStrMat,ZMat,ZpMat,u,zpzMat,keyU,nCols,ycorr,varE,varU,scaleZ,dfZ)
        #for each random effect
        for zSet in keys(zpzMat)
		if isa(zSet,Tuple)
			uPos = keyU[zSet]
			ycorr .+= ZMat[zSet]*u[uPos,1:nCols[zSet]]
			u[uPos,1:nCols[zSet]]  .= sampleU(zSet,iStrMat[zSet],uPos,ZMat[zSet],ZpMat[zSet],zpzMat[zSet],varE,varU[zSet],u[uPos,1:nCols[zSet]],ycorr)
			ycorr .-= ZMat[zSet]*u[uPos,1:nCols[zSet]]
                        varU[zSet] = sampleVarCoVarU(iStrMat[zSet],scaleZ[zSet],dfZ[zSet],u[uPos,1:nCols[zSet]])
		elseif isa(zSet,Symbol) || isa(zSet,Expr)
                	uPos = keyU[zSet]
			ycorr .+= ZMat[zSet]*u[uPos,1:nCols[zSet]]
                	u[uPos,1:nCols[zSet]]  .= sampleU(zSet,iStrMat[zSet],uPos,ZMat[zSet],ZpMat[zSet],zpzMat[zSet],varE,varU[zSet],u[uPos,1:nCols[zSet]],ycorr)
			ycorr .-= ZMat[zSet]*u[uPos,1:nCols[zSet]]		
			varU[zSet] = sampleVarU(iStrMat[zSet],scaleZ[zSet],dfZ[zSet],u[uPos,1:nCols[zSet]])
       		end
		
	 end
end


#sample random marker effects

##### Component-wise, seperated functions for symbol and tuple

function sampleBayesPR!(mSet::Symbol,M,beta,ycorr,varE,varBeta)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize = length(theseLoci)
		lambda = varE/(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = (BLAS.dot(view(M[mSet].data,:,locus),ycorr)) .+ view(M[mSet].rhs,locus)
			lhs = M[mSet].mpm[locus] + lambda
			meanBeta = lhs\rhs
			setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
			BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
		end
@time		@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),regionSize)
	end
end
	
function sampleBayesPR2!(M,beta,ycorr,varE,varBeta)
	for mSet in keys(M)
		local rhs::Float64
		local lhs::Float64
		local meanBeta::Float64
		local lambda::Float64
		for (r,theseLoci) in enumerate(M[mSet].regionArray)
			regionSize = length(theseLoci)
			lambda = varE/(varBeta[mSet][r])
			for locus in theseLoci::UnitRange{Int64}
				BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
				rhs = (BLAS.dot(view(M[mSet].data,:,locus),ycorr)) .+ view(M[mSet].rhs,locus)
				lhs = M[mSet].mpm[locus] + lambda
				meanBeta = lhs\rhs
				setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
				BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			end
	@time		@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),regionSize)
		end
	end
end


##### Component-wise, seperated functions for symbol and tuple
function sampleBayesPR!(mSet::Tuple,MMat,nowMp,beta,mpmMat,betaPos,regionsMat,regions,ycorr,varE,varBeta,scaleMNow,dfMNow,ssRHS)
	for r in 1:regions
		theseLoci = regionsMat[r]
		regionSize = length(theseLoci)
		invB = inv(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			RHS = zeros(size(invB,1))	
			ycorr .+= MMat[locus]*getindex.(beta[betaPos],locus)				
			RHS = ((nowMp[locus]*ycorr)./varE) .+ view(ssRHS,locus) 
			invLHS::Array{Float64,2} = inv((mpmMat[locus]./varE) .+ invB)
			meanBETA::Array{Float64,1} = invLHS*RHS
			setindex!.(beta[betaPos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),locus)
			ycorr .-= MMat[locus]*getindex.(beta[betaPos],locus)	
		end
		varBeta[mSet][r] = sampleVarCovBetaPR(scaleMNow,dfMNow,reduce(hcat,getindex.(beta[betaPos],Ref(theseLoci))),regionSize)
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
	return (scale_ranVar*df_ranVar + effVec'*iMat*effVec)/rand(Chisq(df_ranVar + n))
end

function sampleCoVarU(iMat,scale_ranVar,df_ranVar,effVec)
	n = size(iMat,2)
        return rand(InverseWishart(df_ranVar + n, effVec*iMat*effVec' + scale_ranVar))
end

#sample marker variances
function sampleVarBetaPR(scalem,dfm,whichLoci,regionSize)::Float64
	return (scalem*dfm + BLAS.dot(whichLoci,whichLoci)) / rand(Chisq(dfm + regionSize))
end

function sampleVarCovBetaPR(scalem,dfm,whichLoci,regionSize)
	Sb = whichLoci'whichLoci
	return rand(InverseWishart(dfm + regionSize, scalem + Sb))
end

#Sample residual variance
function sampleVarE(df_e,S_e,yCorVec,nRecords)
    return (df_e*S_e + BLAS.dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords))
end



end
