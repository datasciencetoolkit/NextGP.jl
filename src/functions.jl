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

function sampleBayesPR!(mSet::Symbol,M::Dict,beta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize::Int64 = length(theseLoci)
		lambda = varE/(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rhs = BLAS.dot(view(M[mSet].data,:,locus),ycorr) .+ view(M[mSet].rhs,locus)
			lhs = M[mSet].mpm[locus] + lambda
			meanBeta = lhs\rhs
			setindex!(beta[M[mSet].pos],sampleBeta(meanBeta, lhs, varE),locus)
			BLAS.axpy!(-1.0*getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
		end
		@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],theseLoci),regionSize)
	end
end
	
##### Component-wise, seperated functions for symbol and tuple
function sampleBayesPR!(mSet::Tuple,M::Dict,beta::Vector,ycorr::Vector{Float64},varE::Float64,varBeta::Dict)
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize = length(theseLoci)
		invB = inv(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			RHS = zeros(size(invB,1))	
			ycorr .+= M[mSet].data[locus]*getindex.(beta[M[mSet].pos],locus)
			println("size M[mSet].Mp[locus]: $(M[mSet].Mp[1]./varE), size ycorr: $(ycorr), size rhs: $(view(M[mSet].rhs,1))")
			println("M[mSet].Mp[locus]y: $((M[mSet].Mp[1]*ycorr)./varE), size ycorr: $ycorr, rhs: $(view(M[mSet].rhs,1))")
			RHS = ((M[mSet].Mp[locus]*ycorr)./varE) .+ view(M[mSet].rhs,locus) 
			invLHS::Array{Float64,2} = inv((M[mSet].mpm[locus]./varE) .+ invB)
			meanBETA::Array{Float64,1} = invLHS*RHS
			setindex!.(beta[M[mSet].pos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),locus)
			ycorr .-= M[mSet].data[locus]*getindex.(beta[M[mSet].pos],locus)	
		end
		@inbounds varBeta[mSet][r] = sampleVarCovBetaPR(M[mSet].scale,M[mSet].df,reduce(hcat,getindex.(beta[M[mSet].pos],Ref(theseLoci))),regionSize)
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
