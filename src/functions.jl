module functions

using SpecialFunctions
using Distributions, LinearAlgebra, DataStructures
using StatsBase
#using Printf
#using CSV
#using DataFrames
#using DataStructures

include("outFiles.jl")

export sampleVarE!
export sampleX!
export sampleBayesPR!,sampleBayesB!,sampleBayesC!,sampleBayesR!,sampleBayesRCπ!,sampleBayesLV!,sampleBayesRCplus!
export sampleZ!

#Sampling fixed effects

### NEW, Wang's trick

function sampleb!(xSet::Union{Symbol,Tuple},X::OrderedDict,b::Vector,ycorr::Vector,iVarE)
	#iVarE = inv(varE)
	bVec = deepcopy(b[X[xSet].pos])
	Yi = (X[xSet].Xp*ycorr).*iVarE #computation of X'ycorr*iVarE for ALL  rhsb
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

function sampleb!(xSet::Union{Symbol,Tuple},X::OrderedDict,b::Vector,ycorr::Matrix,iVarE)
	bVec = deepcopy(b[X[xSet].pos])
	for j in 1:X[xSet].nCol
		rhsj0 = []
		Yj = [BLAS.dot(X[xSet].data[j][:,t],sum(ycorr .* iVarE[[t],:],dims=2)) for t in 1:size(ycorr,2)]
        	bVec[j] .= 0.0 #also excludes the effect from iMat! Nice trick.
		foreach(k -> push!(rhsj0,k.*iVarE), X[xSet].XpX[[j],:])
		rr = hcat(rhsj0...)*vec(hcat(bVec...))
		rhsj = Yj - rr
                lhsj = getindex(X[xSet].XpX,j,j).*iVarE
		invLhsj = inv(lhsj)
                meanj = invLhsj*rhsj
                bVec[j] .= reshape(rand(MvNormal(vec(meanj),convert(Array,Symmetric(invLhsj)))),1,size(ycorr,2)) # i convert vector to matrix here
        end
	return bVec
end

# NEW with D and with Wang's Trick (ySet::Symbol)
function sampleX!(xSet::Union{Symbol,Tuple},X::OrderedDict,b::Vector,ycorr::Vector,varE::Dict,ySet::Symbol)
	iVarE = inv(varE[ySet])
	if X[xSet].nCol==1
		ycorr    .+= X[xSet].data .* b[X[xSet].pos]
		rhs      = dot(X[xSet].Xp,ycorr).*iVarE .+ X[xSet].rhs
		lhs      = X[xSet].xpx .*iVarE .+ X[xSet].lhs
		meanMu   = lhs\rhs		
                b[X[xSet].pos] .= rand(Normal(meanMu[],sqrt(inv(lhs[]))))
		ycorr    .-= X[xSet].data .* b[X[xSet].pos]
	else
		ycorr    .+= X[xSet].data*b[X[xSet].pos]
		b[X[xSet].pos] .= sampleb!(xSet,X,b,ycorr,iVarE)
		ycorr    .-= X[xSet].data*b[X[xSet].pos]
	end
end

# NEW with D and with Wang's Trick (ySet::Tuple)
function sampleX!(xSet::Union{Symbol,Tuple},X::OrderedDict,b::Vector,ycorr::Matrix,varE::Dict,ySet::Tuple)
	iVarE = inv(varE[ySet])
	for i in 1:X[xSet].nCol
		ycorr .+= X[xSet].data[i] .* b[X[xSet].pos][i]
	end
	b[X[xSet].pos] .= sampleb!(xSet,X,b,ycorr,iVarE)
	for i in 1:X[xSet].nCol
		ycorr .-= X[xSet].data[i] .* b[X[xSet].pos][i]
	end
end

## Single-trait
#sample random effects
#Single U, D, no-correlation
function sampleU(zSet::Union{Expr,Symbol},Z::OrderedDict,iVarE::Float64,varU::Dict,u::Vector,ycorr::Vector{Float64})
	uVec = deepcopy(u[Z[zSet].pos])
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

#Multiple U correlated, no D
function sampleU(zSet::Tuple,Z::OrderedDict,iVarE::Float64,varU::Dict,u::Vector,ycorr::Vector{Float64})
	uVec = deepcopy(u[Z[zSet].pos])
	nCol = size(uVec,2)
	iVarU = inv(varU[zSet])
	for i in 1:nCol
		setindex!(uVec,[0;0],:,i)
		Yi = Z[zSet].Zp[i]*ycorr.*iVarE	
		rhsU = Yi - kron(view(Z[zSet].iVarStr,[i],:),iVarU)*vec(uVec)
                invLhsU = inv((getindex(Z[zSet].zpz,i).*iVarE) + (view(Z[zSet].iVarStr,i,i).*iVarU))
                meanU = invLhsU*rhsU
		setindex!(uVec,rand(MvNormal(meanU,convert(Array,Symmetric(invLhsU)))),:,i)
        end
	return uVec
end

#Single U Main with D
function sampleZ!(zSet::Union{Expr,Symbol},Z::OrderedDict,u::Vector,ycorr::Vector{Float64},varE::Dict,ySet::Symbol,varU::Dict)
        #for each random effect
	iVarE = inv(varE[ySet])
	ycorr .+= Z[zSet].data*u[Z[zSet].pos]'
        u[Z[zSet].pos] .= sampleU(zSet,Z,iVarE,varU,u,ycorr)
	ycorr .-= Z[zSet].data*u[Z[zSet].pos]'		
	varU[zSet] = sampleVarU(Z[zSet].iVarStr,Z[zSet].scale,Z[zSet].df,u[Z[zSet].pos])		
end

#Multiple U correlated, no D
function sampleZ!(zSet::Tuple,Z::OrderedDict,u::Vector,ycorr::Vector{Float64},varE::Dict,ySet::Symbol,varU::Dict)
	iVarE = inv(varE[ySet])
	nCol = size(u[Z[zSet].pos],2)
	for i in 1:nCol
		ycorr .+= Z[zSet].data[i]*getindex(u[Z[zSet].pos],:,i)
	end
	u[Z[zSet].pos] .= sampleU(zSet,Z,iVarE,varU,u,ycorr)
	varU[zSet] = sampleCoVarU(Z[zSet].iVarStr,Z[zSet].scale,Z[zSet].df,u[Z[zSet].pos])
	for i in 1:nCol
		ycorr .-= Z[zSet].data[i]*getindex(u[Z[zSet].pos],:,i)
	end
end



#sample random marker effects

##### Component-wise, seperated functions for symbol and tuple

function sampleBayesPR!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	iVarE = 1/varE[ySet]
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
		@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale[],M[mSet].df[],getindex(beta[M[mSet].pos],theseLoci),regionSize)
	end
	if (M[mSet].params==true) && (length(varBeta[mSet]) > 1)
		dfIG,scaleIG = sampleScaleDFofVar(varBeta[mSet])
		dfScaleInvChi = 2*dfIG
		scaleScaleInvChi = scaleIG/dfIG
		setindex!(M[mSet].scale, scaleScaleInvChi,1)
		setindex!(M[mSet].df, dfScaleInvChi,1)
	elseif (M[mSet].params==true) && (length(varBeta[mSet]) < 2)
	else nothing #println("What is WRONG???")
	end
end
	
##### Component-wise, seperated functions for symbol and tuple
function sampleBayesPR!(mSet::Tuple,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
	iVarE = 1/varE[ySet]
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize = length(theseLoci)
		invB = inv(varBeta[mSet][r])
		for locus in theseLoci::UnitRange{Int64}
			ycorr .+= M[mSet].data[locus]*getindex(beta[M[mSet].pos],:,locus)
			RHS = ((getindex(M[mSet].Mp,locus)*ycorr).*iVarE)
			invLHS::Array{Float64,2} = inv((getindex(M[mSet].mpm,locus).*iVarE) .+ invB)
			meanBETA::Array{Float64,1} = invLHS*RHS
			beta[M[mSet].pos][:,locus] .= rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS))))
			#setindex!(beta[M[mSet].pos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),locus)
			ycorr .-= M[mSet].data[locus]*getindex(beta[M[mSet].pos],:,locus)	
		end
		@inbounds varBeta[mSet][r] = sampleVarCovBetaPR(M[mSet].scale[],M[mSet].df[],getindex(beta[M[mSet].pos],:,theseLoci),regionSize)
	end
end

##### MULTI-TRAIT BAYESPR###########
function sampleBayesPR!(mSet::Tuple,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Matrix{Float64},varE::Dict,varBeta::Dict,ySet::Tuple)
	iVarE = inv(varE[ySet])
	for (r,theseLoci) in enumerate(M[mSet].regionArray)
		regionSize = length(theseLoci)
		invB = inv(varBeta[mSet][r]) 
		for locus in theseLoci::UnitRange{Int64}
			if isa(mSet,Tuple{Vararg{Symbol}}) #could also be Tuple{Vararg{Tuple{Vararg{Symbol}}}} which i will adapt later
				#for i in 1:length(ySet)
				#	ycorr[:,i] .+= M[mSet].data[locus][:,i]*getindex(beta[M[mSet].pos],i,locus)
				#end
				ycorr .+= M[mSet].data[locus].*beta[M[mSet].pos][:,locus]'
			end
			RHS = sum(((getindex(M[mSet].Mp,locus)*ycorr).*iVarE),dims=2)
			invLHS::Array{Float64,2} = inv((getindex(M[mSet].mpm,locus).*iVarE) .+ invB)
			meanBETA = vec(invLHS*RHS)			
			setindex!(beta[M[mSet].pos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),:,locus)
			if isa(mSet,Tuple{Vararg{Symbol}})
				#for i in 1:length(ySet)
				#	ycorr[:,i] .-= M[mSet].data[locus][:,i]*getindex(beta[M[mSet].pos],i,locus)
				#end
				ycorr .-= M[mSet].data[locus].*beta[M[mSet].pos][:,locus]'
			end
		end
		#@inbounds varBeta[mSet][r] = M[mSet].params==true ? sampleVarCovBetaPR(M[mSet].scale[],M[mSet].df[],getindex(beta[M[mSet].pos],:,theseLoci),regionSize) : sampleVarCovBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],:,theseLoci),regionSize) 
		@inbounds varBeta[mSet][r] = M[mSet].params == true ? sampleHierarchicalSIW!(M[mSet], delta, mSet, getindex(beta[M[mSet].pos], :, theseLoci), regionSize) : sampleVarCovBetaPR(M[mSet].scale, M[mSet].df, getindex(beta[M[mSet].pos], :, theseLoci), regionSize)
	end
	if (M[mSet].params==true) #&& (length(varBeta[mSet]) > 1) #estimate in all cases
		#dfIG,scaleIG = sampleScaleDFofVar(varBeta[mSet])
		#dfScaleInvChi = 2*dfIG
		#scaleScaleInvChi = scaleIG/dfIG
		#setindex!(M[mSet].scale, scaleScaleInvChi, 1)
		#setindex!(M[mSet].df, dfScaleInvChi, 1)
		#scaleEst = rand(InverseWishart(5, mean(inv.(varBeta[mSet])) + inv(M[mSet].priorScale)))
	elseif (M[mSet].params==true) && (length(varBeta[mSet]) < 2)
	else nothing #println("What is WRONG???")
	end
end


function sampleBayesB!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	nLoci = 0
	iVarE = 1/varE[ySet]
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		iVarBeta = 1/varBeta[mSet][r]
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rrr = BLAS.dot(view(M[mSet].data,:,locus),ycorr) #not rhs!
			v0 = getindex(M[mSet].mpm,locus)*varE[ySet]
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
				@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale[],M[mSet].df[],getindex(beta[M[mSet].pos],theseLoci),1)
			else 
				setindex!(beta[M[mSet].pos],0.0,locus)
				setindex!(delta[M[mSet].pos],0,locus)
				@inbounds varBeta[mSet][r] = sampleVarBetaPR(M[mSet].scale[],M[mSet].df[],[0.0],0)
			end
		end
	end
	if M[mSet].estPi == true 
		piIn = samplePi(nLoci,M[mSet].dims[2]) #probability of in
		M[mSet].piHat .= [1.0-piIn piIn]
		M[mSet].logPi .= log.([1.0-piIn piIn])
	end
	if (M[mSet].params==true) && (length(varBeta[mSet]) > 1)
		inMarkers = findall(>(0),delta[M[mSet].pos])
		inMarkers = getindex.(inMarkers,2)
		dfIG,scaleIG = sampleScaleDFofVar(getindex(varBeta[mSet],inMarkers))
		#dfIG,scaleIG = sampleScaleDFofVar(varBeta[mSet])
		dfScaleInvChi = 2*dfIG
		scaleScaleInvChi = scaleIG/dfIG
		setindex!(M[mSet].scale, scaleScaleInvChi,1)
		setindex!(M[mSet].df, dfScaleInvChi,1)
	elseif (M[mSet].params==true) && (length(varBeta[mSet]) < 2)
	else nothing #println("What is WRONG???")
	end
end


##### MULTI-TRAIT BAYESB###########
function sampleBayesB!(mSet::Tuple,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Matrix{Float64},varE::Dict,varBeta::Dict,ySet::Tuple)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	storeCount = zeros(length(M[mSet].gammaComb))
	iVarE = inv(varE[ySet])
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		for locus in theseLoci::UnitRange{Int64}

			iVarBeta = inv(varBeta[mSet][locus])
			
			if isa(mSet,Tuple{Vararg{Symbol}}) #could also be Tuple{Vararg{Tuple{Vararg{Symbol}}}} which i will adapt later
				#for i in 1:length(ySet)
				#	ycorr[:,i] .+= M[mSet].data[locus][:,i]*M[mSet][:gammaHat][locus][i]*getindex(beta[M[mSet].pos],i,locus)
				#end
				ycorr .+= M[mSet].data[locus].*M[mSet][:gammaHat][locus].*beta[M[mSet].pos][:,locus]'

			end

			tempGammaProb = zeros(length(M[mSet].gammaComb)) #to store prob
			
			for thisGamma in 1:length(M[mSet].gammaComb)
               			C = M[mSet].deltaComb[thisGamma]*getindex(M[mSet].mpm,locus)*M[mSet].deltaComb[thisGamma].*iVarE .+ iVarBeta
                		invC = inv(C)
				rhsReg = vec(sum(((M[mSet][:deltaHat][locus]*getindex(M[mSet].Mp,locus)*ycorr).*iVarE),dims=2))
				logL = -(0.5)*log(det(invC) + (rhsReg'*invC*rhsReg)) + log(M[mSet].piHat[thisGamma])
				tempGammaProb[thisGamma] = exp(logL)
            		end
			prob4Region = tempGammaProb./sum(tempGammaProb)
			myDelta = rand(Categorical(prob4Region))
			
            		M[mSet][:gammaHat][locus] = M[mSet].gammaComb[myDelta]
            		M[mSet][:deltaHat][locus] = M[mSet].deltaComb[myDelta]
            		storeCount[myDelta] += 1.0

			RHS = vec(sum(((M[mSet][:deltaHat][locus]*getindex(M[mSet].Mp,locus)*ycorr).*iVarE),dims=2))
			invLHS::Array{Float64,2} = inv(M[mSet][:deltaHat][locus]*(getindex(M[mSet].mpm,locus)*M[mSet][:deltaHat][locus].*iVarE) .+ iVarBeta)
			meanBETA = vec(invLHS*RHS)			
			setindex!(beta[M[mSet].pos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),:,locus)
			setindex!(delta[M[mSet].pos],M[mSet][:gammaHat][locus],:,locus) ###DOUBLE WORK

			if isa(mSet,Tuple{Vararg{Symbol}}) #could also be Tuple{Vararg{Tuple{Vararg{Symbol}}}} which i will adapt later
				#for i in 1:length(ySet)
				#	ycorr[:,i] .-= M[mSet].data[locus][:,i]*M[mSet][:gammaHat][locus][i]*getindex(beta[M[mSet].pos],i,locus)
				#end
				ycorr .-= M[mSet].data[locus].*M[mSet][:gammaHat][locus].*beta[M[mSet].pos][:,locus]'
			end
		end
		@inbounds varBeta[mSet][r] = M[mSet].params==true ? sampleVarCovBetaPR(M[mSet].scale[],M[mSet].df[],getindex(beta[M[mSet].pos],:,theseLoci),1) : sampleVarCovBetaPR(M[mSet].scale,M[mSet].df,getindex(beta[M[mSet].pos],:,theseLoci),1) 
	end
	#println("storeCount: $storeCount")
	if M[mSet].estPi == true 
		M[mSet].piHat .= rand(Dirichlet(storeCount.+1))
	end
	if (M[mSet].params==true) #&& (length(varBeta[mSet]) > 1) #estimate in all cases
		#println("$(M[mSet].scale[]),$(M[mSet].df[])")
		#inMarkers = findall([any(c.!=0) for c in eachcol(delta[M[mSet].pos])])
		inMarkers = findall([any(c.!=0) for c in M[mSet][:gammaHat]])
		#println("N inmarkers = $(length(inMarkers))")
		scaleEst = rand(InverseWishart(5, convert(Array,Symmetric(inv(mean(getindex(varBeta[mSet],inMarkers))) + inv(M[mSet].priorScale)))))
		#NEW
		#println("new scale: $(scaleEst)")
		setindex!(M[mSet].scale, scaleEst,1)

		tempdf = 2 #one for each for the 4+1 
		for m in 1:length(mSet)
			inMarkers = findall(>(0),getindex(beta[M[mSet].pos],m,:))
			dfIG,scaleIG = sampleScaleDFofVar(getindex.(getindex(varBeta[mSet],inMarkers),m,m))
			dfScaleInvChi = 2*dfIG
			tempdf += dfScaleInvChi
		end
		#println("new DF: $(tempdf/length(mSet))")
		setindex!(M[mSet].df,tempdf/length(mSet),1)

	elseif (M[mSet].params==true) && (length(varBeta[mSet]) < 2)
		nothing
	else 
		#println("$(M[mSet].scale),$(M[mSet].df)")
		nothing
	end
end


function sampleBayesC!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	nLoci = 0
	iVarE = 1/varE[ySet]
	iVarBeta = 1/varBeta[mSet][1]
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesC, so r=locus
		for locus in theseLoci::UnitRange{Int64}
			BLAS.axpy!(getindex(beta[M[mSet].pos],locus),view(M[mSet].data,:,locus),ycorr)
			rrr = BLAS.dot(view(M[mSet].data,:,locus),ycorr)
			v0 = getindex(M[mSet].mpm,locus)*varE[ySet]
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
	@inbounds varBeta[mSet][1] = sampleVarBetaPR(M[mSet].scale[],M[mSet].df[],beta[M[mSet].pos],nLoci)
	if M[mSet].estPi == true 
		piIn = samplePi(nLoci,M[mSet].dims[2]) #probability of in
		M[mSet].piHat .= [1.0-piIn piIn]
		M[mSet].logPi .= log.([1.0-piIn piIn])
	end
	#println("scale before: $(M[mSet].scale)")
	#println(M[mSet].df," ", varBeta[mSet]," ", length(M[mSet].mpm))
	M[mSet].params==true ? setindex!(M[mSet].scale, sampleScaleOfVar(M[mSet].df[],varBeta[mSet],M[mSet].nVarCov), 1) : nothing
	#println("scale after: $(M[mSet].scale)")
end


##### MULTI-TRAIT BAYESC###########
function sampleBayesC!(mSet::Tuple,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Matrix{Float64},varE::Dict,varBeta::Dict,ySet::Tuple)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64
	storeCount = zeros(length(M[mSet].gammaComb))
	iVarE = inv(varE[ySet])
	for (r,theseLoci) in enumerate(M[mSet].regionArray) #theseLoci is always as 1:1,2:2 for BayesB, so r=locus
		for locus in theseLoci::UnitRange{Int64}

			iVarBeta = inv(varBeta[mSet][1])
			
			if isa(mSet,Tuple{Vararg{Symbol}}) #could also be Tuple{Vararg{Tuple{Vararg{Symbol}}}} which i will adapt later
				#for i in 1:length(ySet)
				#	ycorr[:,i] .+= M[mSet].data[locus][:,i]*M[mSet][:gammaHat][locus][i]*getindex(beta[M[mSet].pos],i,locus)
				#end
				ycorr .+= M[mSet].data[locus].*M[mSet][:gammaHat][locus].*beta[M[mSet].pos][:,locus]'
			end

			tempGammaProb = zeros(length(M[mSet].gammaComb)) #to store prob
			
			for thisGamma in 1:length(M[mSet].gammaComb)
               			C = M[mSet].deltaComb[thisGamma]*getindex(M[mSet].mpm,locus)*M[mSet].deltaComb[thisGamma].*iVarE .+ iVarBeta
                		invC = inv(C)
				rhsReg = vec(sum(((M[mSet][:deltaHat][locus]*getindex(M[mSet].Mp,locus)*ycorr).*iVarE),dims=2))
				logL = -(0.5)*log(det(invC) + (rhsReg'*invC*rhsReg)) + log(M[mSet].piHat[thisGamma])
				tempGammaProb[thisGamma] = exp(logL)
            		end
			prob4Region = tempGammaProb./sum(tempGammaProb)
			myDelta = rand(Categorical(prob4Region))
			
            		M[mSet][:gammaHat][locus] = M[mSet].gammaComb[myDelta] #vector of diagonals
            		M[mSet][:deltaHat][locus] = M[mSet].deltaComb[myDelta] #matrix form
            		storeCount[myDelta] += 1.0

			RHS = vec(sum(((M[mSet][:deltaHat][locus]*getindex(M[mSet].Mp,locus)*ycorr).*iVarE),dims=2))
			invLHS::Array{Float64,2} = inv(M[mSet][:deltaHat][locus]*(getindex(M[mSet].mpm,locus)*M[mSet][:deltaHat][locus].*iVarE) .+ iVarBeta)
			meanBETA = vec(invLHS*RHS)			
			setindex!(beta[M[mSet].pos],rand(MvNormal(meanBETA,convert(Array,Symmetric(invLHS)))),:,locus)
			setindex!(delta[M[mSet].pos],M[mSet][:gammaHat][locus],:,locus) ###DOUBLE WORK

			if isa(mSet,Tuple{Vararg{Symbol}}) #could also be Tuple{Vararg{Tuple{Vararg{Symbol}}}} which i will adapt later
				#for i in 1:length(ySet)
				#	ycorr[:,i] .-= M[mSet].data[locus][:,i]*M[mSet][:gammaHat][locus][i]*getindex(beta[M[mSet].pos],i,locus)
				#end
				ycorr .-= M[mSet].data[locus].*M[mSet][:gammaHat][locus].*beta[M[mSet].pos][:,locus]'
			end
		end
	end
	@inbounds varBeta[mSet][1] = M[mSet].params==true ? sampleVarCovBetaPR(M[mSet].scale[],M[mSet].df[],beta[M[mSet].pos],M[mSet].dims[2]) : sampleVarCovBetaPR(M[mSet].scale,M[mSet].df,beta[M[mSet].pos],M[mSet].dims[2]) 

	#println("storeCount: $storeCount")
	if M[mSet].estPi == true 
		M[mSet].piHat .= rand(Dirichlet(storeCount.+1))
	end
	if (M[mSet].params==true) #&& (length(varBeta[mSet]) > 1) #estimate in all cases
		#println("$(M[mSet].scale[]),$(M[mSet].df[])")
		#inMarkers = findall([any(c.!=0) for c in eachcol(delta[M[mSet].pos])])
		inMarkers = findall([any(c.!=0) for c in M[mSet][:gammaHat]])
		#println("N inmarkers = $(length(inMarkers))")
		scaleEst = rand(InverseWishart(5, convert(Array,Symmetric(inv(varBeta[mSet][1]) + inv(M[mSet].priorScale)))))
		#NEW
		setindex!(M[mSet].scale, scaleEst,1)
		setindex!(M[mSet].df,rand(Uniform(4,5)),1)

	elseif (M[mSet].params==true) && (length(varBeta[mSet]) < 2)
		nothing
	else 
		#println("$(M[mSet].scale),$(M[mSet].df)")
		nothing
	end
end


function sampleBayesR!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
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
	@inbounds varBeta[mSet][1] = sampleVarBetaR(M[mSet].scale[],M[mSet].df,sumS,nNonZero)
	##
	
	if M[mSet].estPi == true 
		piHat = samplePi(nLoci)
		M[mSet].piHat .= piHat
		M[mSet].logPi .= log.(piHat)
	end
	#println("scale before: $(M[mSet].scale)")
	#println(M[mSet].df," ", varBeta[mSet]," ", length(M[mSet].mpm))
	M[mSet].params==true ? setindex!(M[mSet].scale, sampleScaleOfVar(M[mSet].df,varBeta[mSet],M[mSet].nVarCov), 1) : nothing
	#println("scale after: $(M[mSet].scale)")
end

function sampleBayesRCπ!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
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
		@inbounds varBeta[mSet][a] = sampleVarBetaR(M[mSet].scale[],M[mSet].df,sumS[a],nNonZero[a])
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
	#println("scale before: $(M[mSet].scale)")
	#println(M[mSet].df," ", varBeta[mSet]," ", length(M[mSet].mpm))
	M[mSet].params==true ? setindex!(M[mSet].scale, sampleScaleOfVar(M[mSet].df,varBeta[mSet],M[mSet].nVarCov), 1) : nothing
	#println("scale after: $(M[mSet].scale)")
end

function sampleBayesRCplus!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
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
		@inbounds varBeta[mSet][a] = sampleVarBetaR(M[mSet].scale[],M[mSet].df,sumS[a],nNonZero[a])
	end
	
	#Estimate both Prob (p) and pi π
	if M[mSet].estPi == true
		for a in 1:nAnnot
			piHat = samplePi(nLoci[a,:])
			M[mSet].piHat[a] = piHat
			M[mSet].logPi[a] = log.(piHat)
		end
	end
	#println("scale before: $(M[mSet].scale)")
	#println(M[mSet].df," ", varBeta[mSet]," ", length(M[mSet].mpm))
	M[mSet].params==true ? setindex!(M[mSet].scale, sampleScaleOfVar(M[mSet].df,varBeta[mSet],M[mSet].nVarCov), 1) : nothing
	#println("scale after: $(M[mSet].scale)")
end

function sampleBayesLV!(mSet::Symbol,M::OrderedDict,beta::Vector,delta::Vector,ycorr::Vector{Float64},varE::Dict,varBeta::Dict,ySet::Symbol)
	local rhs::Float64
	local lhs::Float64
	local meanBeta::Float64
	local lambda::Float64

	var_var = M[mSet].varZeta[]
	
	nLoci = 0
	iVarE = 1/varE[ySet]
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
	#println("trapped: $(trapped/(trapped+notTrapped))")

	rhsC = M[mSet].covariatesT*M[mSet].logVar
	meanC = M[mSet].iCpC*rhsC
	M[mSet].c .= rand(MvNormal(vec(meanC),convert(Array,Symmetric(M[mSet].iCpC*var_var))))
	M[mSet].SNPVARRESID .= M[mSet].logVar .- M[mSet].covariates*M[mSet].c

	#sample varZeta
	if isa(M[mSet].estVarZeta,Float64)
		setindex!(M[mSet].varZeta,M[mSet].estVarZeta*var(M[mSet].logVar),1)
	elseif M[mSet].estVarZeta == false
		nothing		
	elseif M[mSet].estVarZeta == true
		setindex!(M[mSet].varZeta,var(M[mSet].SNPVARRESID),1)
	end	
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
	Sb = whichLoci*whichLoci'
	return rand(InverseWishart(dfm + regionSize, scalem + Sb))
end
				
function sampleVarBetaR(scalem,dfm,sumS,nLoci)::Float64
	return (scalem*dfm + sumS) / rand(Chisq(dfm + nLoci))
end

#Sample residual variance
function sampleVarE!(eSet::Symbol,E,varE,yCorVec,nRecords)
	#println("eSet in varE $eSet")
	#println("E in varE $E")
	if E[eSet].str == "D"
		varE[eSet] = (E[eSet].df*E[eSet].scale + sum(E[eSet].iVarStr.*(yCorVec.^2)))/rand(Chisq(E[eSet].df + nRecords))
	else varE[eSet] = (E[eSet].df*E[eSet].scale + BLAS.dot(yCorVec,yCorVec))/rand(Chisq(E[eSet].df + nRecords))
	end
end

#Sample residual variance
function sampleVarE!(eSet::Tuple,E,varE,yCorVec,nRecords)
	Se = yCorVec'yCorVec
	varE[eSet] = rand(InverseWishart(E[eSet].df + nRecords, convert(Array,Symmetric(E[eSet].scale + Se))))
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

function sampleScaleOfVar(df,var,k)
	#length(var) should be nVarComp
	#+1 is for Gamma(1,1) prior
	#println("df: $df 1/var: $(1.0./var) lengthVar: $(length(var))")
	return rand(Gamma(0.5*df*length(var)+1,1/(1+0.5*df*sum(1.0./var))))
end

function sampleScaleDFofVar(varVec;abcde=[4.0,0.01,0.01,0.01,0.01],crit=10e-6)
	a,b,c,d,e = abcde
	n = length(varVec)
	mu = mean(varVec)
	v = var(varVec,corrected=true)
	alpha = ((mu^2)/v)+2
	eHat = e+sum(1.0 ./ varVec)
	logAHat = log(a) + sum(log.(varVec))
	bHat = b + n
	cHat = c + n
	C = log(eHat)
	iter=0	
	while true
		iter+=1
		#println(iter)
		global alphaHat = invdigamma((-logAHat+cHat*(log(d+n*alpha)-C))/(bHat))
		if abs(alphaHat - alpha) >= crit
			alpha = alphaHat
		else
			break
		end
	end
	dHat = d + n*alphaHat
	betaHat = dHat/eHat
	return alphaHat,betaHat
end



using LinearAlgebra
using Distributions
using SpecialFunctions: loggamma

"""
    sampleHierarchicalSIW!(modelRegion, delta_chain, mSet, B_region, regionSize)

A single-line replacement wrapper that runs the fully hierarchical Scaled Inverse Wishart 
updates for a specific genomic region and returns the reconstructed covariance matrix.
"""
function sampleHierarchicalSIW!(modelRegion, delta_chain::Vector, mSet::Tuple, B_region::Matrix{Float64}, regionSize::Int64)
    t = size(B_region, 1)
    
    # 1. Safely extract mutable fields (handles both raw values or Ref containers)
    S_g  = modelRegion.scale isa Ref ? modelRegion.scale[] : modelRegion.scale
    nu_g = modelRegion.df isa Ref    ? modelRegion.df[]    : modelRegion.df
    
    # 2. Extract or define hyper-priors/tuning values dynamically
    nu0_W     = hasproperty(modelRegion, :nu0_W)     ? modelRegion.nu0_W     : t + 2.0
    S0_W_inv  = hasproperty(modelRegion, :S0_W_inv)  ? modelRegion.S0_W_inv  : inv(Matrix(I(t)) * (t + 2.0))
    lambda_nu = hasproperty(modelRegion, :lambda_nu) ? modelRegion.lambda_nu : 0.2
    prop_sd   = hasproperty(modelRegion, :prop_sd)   ? modelRegion.prop_sd   : 0.2
    mu_0      = 0.0
    sigma0_2  = 100.0
    alpha_0   = 4.0
    beta_0    = 4.0

    # 3. Handle trait log-normal scales and structural vectors
    delta_g = delta_chain[modelRegion.pos] 
    m_g     = hasproperty(modelRegion, :m_g)  ? modelRegion.m_g  : zeros(t)
    s2_g    = hasproperty(modelRegion, :s2_g) ? modelRegion.s2_g : ones(t) * 2.0
    
    # 4. Step A: Update trait scale standard deviations (delta_g) via MH step
    Qg_old = hasproperty(modelRegion, :Qg) ? modelRegion.Qg : inv(S_g * (nu_g - t - 1))
    _siw_update_scales!(delta_g, B_region, Qg_old, m_g, s2_g, prop_sd)
    
    # 5. Step B: Strip out scales to generate the unscaled scatter matrix and sample Qg
    inv_delta = 1.0 ./ delta_g
    B_unscaled = B_region .* inv_delta
    scatter_unscaled = B_unscaled * B_unscaled'
    
    Scale_InvWish = Symmetric(S_g + scatter_unscaled)
    L_Q = cholesky(inv(Scale_InvWish)).L
    nu_Q = nu_g + regionSize
    A_Q = zeros(t, t)
    for i in 1:t
        A_Q[i, i] = sqrt(rand(Chi(nu_Q - i + 1)))
        for j in 1:(i-1); A_Q[i, j] = randn(); end
    end
    Qg_new = inv(Symmetric((L_Q * A_Q) * (L_Q * A_Q)'))
    if hasproperty(modelRegion, :Qg); modelRegion.Qg .= Qg_new; end
    
    # 6. Step C: Update Wishart scale matrix (S_g) and Shifted Exponential df (nu_g)
    nu_g_new, _ = _siw_update_structure!(Qg_new, nu_g, S_g, nu0_W, S0_W_inv, lambda_nu, prop_sd)
    
    # 7. Step D: Update log-normal hyperparameters (m_g and s2_g) online
    for k in 1:t
        log_d = log(delta_g[k])
        V_m = 1.0 / (1.0 / sigma0_2 + 1.0 / s2_g[k])
        m_g[k] = V_m * (mu_0 / sigma0_2 + log_d / s2_g[k]) + sqrt(V_m) * randn()
        s2_g[k] = rand(InverseGamma(alpha_0 + 0.5, beta_0 + 0.5 * (log_d - m_g[k])^2))
    end
    
    # 8. Write modified scalars back to your reference container variables
    if modelRegion.df isa Ref; modelRegion.df[] = nu_g_new; else modelRegion.df = nu_g_new; end
    # S_g was mutated explicitly in place inside _siw_update_structure!

    # 9. Return the newly generated covariance matrix directly to varBeta[mSet][r]
    return Diagonal(delta_g) * Qg_new * Diagonal(delta_g)
end

# Internal Subroutine: Updates scale standard deviation vector (delta_g)
function _siw_update_scales!(delta_g::Vector{Float64}, B::Matrix{Float64}, Qg::Matrix{Float64}, m_g::Vector{Float64}, s2_g::Vector{Float64}, prop_sd::Float64)
    t, p = size(B)
    BtB = B * B'
    
    function log_target_delta(dg_vec::Vector{Float64})
        D = Diagonal(dg_vec)
        G = D * Qg * D
        G_fact = try cholesky(Symmetric(G)) catch; return -Inf end
        
        log_lik = -0.5 * p * logdet(G_fact) - 0.5 * tr(inv(G_fact) * BtB)
        log_prior = 0.0
        for k in 1:t
            log_prior += -0.5 * (log(dg_vec[k]) - m_g[k])^2 / s2_g[k] - log(dg_vec[k])
        end
        return log_lik + log_prior
    end

    for k in 1:t
        log_delta_curr = log(delta_g[k])
        log_delta_prop = log_delta_curr + prop_sd * randn()
        delta_prop_val = exp(log_delta_prop)
        
        delta_vec_prop = copy(delta_g)
        delta_vec_prop[k] = delta_prop_val
        
        log_alpha = log_target_delta(delta_vec_prop) - log_target_delta(delta_g)
        if log(rand()) < log_alpha
            delta_g[k] = delta_prop_val
        end
    end
end

# Internal Subroutine: Updates scale matrix (S_g) and degrees of freedom (nu_g)
function _siw_update_structure!(Qg::Matrix{Float64}, nu_g::Float64, S_g::Matrix{Float64}, nu0_W::Float64, S0_W_inv::Matrix{Float64}, lambda_nu::Float64, prop_sd::Float64)
    t = size(Qg, 1)
    c_threshold = t - 1.0
    
    # Conjugate Wishart update for S_g matrix elements
    Qg_inv = inv(Symmetric(Qg))
    S_n_inv = S0_W_inv + Qg_inv
    nu_n = nu0_W + nu_g
    L_S = cholesky(inv(Symmetric(S_n_inv))).L
    
    A = zeros(t, t)
    for i in 1:t
        A[i, i] = sqrt(rand(Chi(nu_n - i + 1)))
        for j in 1:(i-1); A[i, j] = randn(); end
    end
    X = L_S * A
    S_g .= X * X' # Mutates elements in place seamlessly

    # Random-Walk Metropolis-Hastings step with Shifted Exponential prior for nu_g
    function log_target_nu(nu::Float64)
        log_prior = log(lambda_nu) - lambda_nu * (nu - c_threshold)
        log_mv_gamma = (t * (t - 1) / 4) * log(pi) + sum(loggamma((nu + 1.0 - i) / 2.0) for i in 1:t)
        log_det_block = (nu / 2.0) * logdet(Symmetric(S_g)) - (nu * t / 2.0) * log(2.0) - log_mv_gamma - (nu / 2.0) * logdet(Symmetric(Qg))
        return log_det_block + log_prior
    end

    nu_prop = nu_g + prop_sd * randn()
    accepted = false
    if nu_prop > c_threshold
        log_alpha = log_target_nu(nu_prop) - log_target_nu(nu_g)
        if log(rand()) < log_alpha
            nu_g = nu_prop
            accepted = true
        end
    end
    return nu_g, accepted
end



							
end
