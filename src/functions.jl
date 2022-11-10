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
function sampleX!(xSet,xMat::Array{Float64, 1},b,ixpx,pos,ycorr,varE,ssRHS)
        #block for each effect
		ycorr    .+= xMat.*b[pos]
		rhs      = xMat'*ycorr .+ ssRHS
		meanMu   = ixpx*rhs			
                b[pos] .= rand(Normal(meanMu[],sqrt((ixpx*varE))[]))
		ycorr    .-= xMat.*b[pos]        
end

function sampleX!(xSet,xMat::Array{Float64, 2},b,ixpx,pos,ycorr,varE,ssRHS)
        #block for each effect
		ycorr    .+= xMat*b[pos]
                rhs      = xMat'*ycorr .+ ssRHS
                meanMu   = ixpx*rhs
		b[pos] .= rand(MvNormal(vec(meanMu),convert(Array,Symmetric(ixpx*varE))))
		ycorr    .-= xMat*b[pos]
end

#sample random effects
function sampleU(zSet::Symbol,iMat,pos,ZComp,ZpComp,zpzComp,varE,varUComp,uVector,ycorr)
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


function sampleBayesPR!(mSet::Symbol,MMat,beta,ycorr,varE,varBeta)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	for r in 1:MMat.nRegions
		theseLoci = MMat.regionsArray[r]
		regionSize = length(theseLoci)
		lambda = varE/(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[MMat.pos],locus),view(MMat.data,:,locus),ycorr)
			rhs = (BLAS.dot(view(MMat.data,:,locus),ycorr)) .+ view(MMat.rhs,locus)
			lhs = MMat.mpm[locus] + lambda
			meanBeta = lhs\rhs
			setindex!(beta[MMat.pos],sampleBeta(meanBeta, lhs, varE),locus)
			BLAS.axpy!(-1.0*getindex(beta[MMat.pos],locus),view(MMat.data,:,locus),ycorr)
		end
		varBeta[mSet][r] = sampleVarBetaPR(MMat.scale,MMat.df,getindex(beta[MMat.pos],theseLoci),regionSize)
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
