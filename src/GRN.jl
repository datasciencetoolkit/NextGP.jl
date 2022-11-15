module GRN

using StatsBase
using LinearAlgebra
using Distributions
using CSV
using DelimitedFiles

export estGRN_MHGibbs

#===
No variable selection
===#
function estGRN_MHGibbs(X,Y,nGenes,SNP4Gene,chainLength,burnIn,outputFreq;startλ1=ones((nGenes^2)-nGenes)*0.0,meanλ1 = 0.0,startΛ2=zeros(nGenes,nGenes*SNP4Gene),priorRes=1.0,outFolder="outMCMC")
    
    #create folder
    folderHandler(outFolder)

    #store results
    these2Keep  = collect((burnIn+outputFreq):outputFreq:chainLength)
    
    #center X. Y corrected later
    Xc = X .- mean(X,dims=2) #SNPs,ind
    nInd  = size(X,2)
    nSNPs = size(X,1) #SNPs,ind
    println("nSNPs = $nSNPs")
    
    nRecords = nGenes*nInd #??????
    
#    nRegions = 1
    
    SNPList = []
    counter = 1
    for i in 1:nGenes
        push!(SNPList,counter:i*SNP4Gene)
        counter = i*SNP4Gene+1
    end
    SNPList = collect.(SNPList)
        
    #initial values
    #Y is genes ind
    μ = vec(mean(Y,dims=2)) # zeros(nGenes)
    Λ2 = startΛ2

    #lambda1 settings
    #starting value, but also current values I  have zeros
    currentλ1 = startλ1
    #to store current Lambda1... diagonal set to zero???
    currentΛ1 = Matrix(I*0.0,nGenes,nGenes)
    #code to get positions of off-diagonals of Lambda1
    get_offdiagPos(A) = [CartesianIndex(ι[2], ι[1]) for ι in CartesianIndices(A) if ι[1] ≠ ι[2]]
    posΛ1 = get_offdiagPos(currentΛ1)
    #Now current Lambda1 is consistent with lambda1
    currentΛ1[posΛ1] .= currentλ1
    println("Starting values of Λ1: $(hcat(currentΛ1...))")
    #n lambda1 coef
    nλ1Coeff = length(posΛ1)
    #allocate memory for diagonal matrix (what about using  diagind()??)
    sum_MpM_I = Matrix(I*1.0,nλ1Coeff,nλ1Coeff)
    LambdaStar_I = Matrix(I*1.0,nGenes,nGenes)
    #acceptance counter
    accept = [0]
    #to store proposed lambda1... diagonal set to zero???
    proposedΛ1 = Matrix(I*0.0,nGenes,nGenes)
    #lambda1 settings ended
    
    #priors
    dfL1Var     = 4.0
    dfEffectVar = 4.0
    dfRes       = 4.0
    means4λ1    = meanλ1
    means4Λ2    = zeros(nGenes) ###should be changed later
    
    #precomputation of vsB for convenience
    ###############################These are arbitrary values.
    varLambda1   = 0.0005            
    varBeta      = fill(0.0005, nGenes) #each gene has its own variance
    varResidual  = priorRes
    
    scaleL1Var       = varLambda1*(dfL1Var-2.0)/dfL1Var
    νS_L1            = scaleL1Var*dfL1Var
    df_L1            = dfL1Var
    
    scaleVar        = varBeta[1]*(dfEffectVar-2.0)/dfEffectVar
    νS_β            = scaleVar*dfEffectVar
    df_β            = dfEffectVar
                
    scaleRes        = varResidual*(dfRes-2.0)/dfRes
    νS_e            = scaleRes*dfRes
    df_e            = dfRes    
    
    #Start analyses

#    yCorr = Y-μ*ones(N)'-Λ2*X
    yCorr = Y.-μ-currentΛ1*Y-Λ2*Xc
                
    #Lambda1 settings
    BIGM = lambda1BIGM(nGenes,yCorr,nInd)
    #pre-compute MpM            
    sum_MpM = BIGM'BIGM

    #Sampling!!!!!!
    for iter in 1:chainLength
                
        #1) sample residual variance
        varE = sampleVarE(νS_e,yCorr,df_e,nRecords) #dot product in the function already vectorizes yCorr.
        
        #2) sample means
        yCorr .+= μ
                    
        rhs      = sum(yCorr,dims=2)
        invLhs   = 1.0/nInd
        meanMu   = rhs.*invLhs
        for g in 1:nGenes
            μ[g] = rand(Normal(meanMu[g],sqrt(invLhs*varE)))
        end
#        μ .= rand(MvNormal(vec(meanMu),invLhs*varE))         
        
        yCorr .-= μ

        #3) sample Λ1. MH
        sampleΛ1!(yCorr,nInd,sum_MpM,BIGM,posΛ1,nλ1Coeff,sum_MpM_I,LambdaStar_I,currentλ1,currentΛ1,accept,proposedΛ1;λ1_0=meanλ1,σ2ϵ=varE,τ2=varLambda1)             
                    
        #4) sample tao1^2
        varLambda1 = sampleVarLambda1(νS_L1,currentλ1.-meanλ1,df_L1,nλ1Coeff)
                    
        #5) sample Λ2. single-site gibbs
        sampleΛ2!(Λ2,Xc,yCorr,varBeta,varE,means4Λ2)
        
        #6) sample SNP variances tao2^2. Could be combined with above function
        for g in 1:nGenes
            varBeta[g] = sampleVarBeta(νS_β,Λ2[g,:],df_β,nSNPs)
        end
                                        
        if iter in these2Keep
            outMCMC(outFolder,"Lambda1",hcat(currentΛ1...))
            outMCMC(outFolder,"varLambda1",hcat(varLambda1...))
            outMCMC(outFolder,"Lambda2",hcat(Λ2...))
            outMCMC(outFolder,"varBeta",hcat(varBeta...))
            outMCMC(outFolder,"varE",varE)
            outMCMC(outFolder,"means",hcat(μ...))
        end
    end
    return accept[]
end


##if really want to use BLAS.axpy!, then yCorr must be made array of arrays
##yCorr = [yCorr[:,1],yCorr[:,2],.......], than BLAS.axpy! will work. Otherwise it wont, as I've experienced before!
function sampleΛ2!(Λ2,Xc,yCorr,σ2τ,σ2ϵ,pMeans)
    nGenes,nSNPs = size(yCorr,1),size(Xc,1) #Xc SNPs,ind;w genes,ind
    for g in 1:nGenes
        α = σ2ϵ/σ2τ[g]
        for q in 1:nSNPs
            yCorr[g,:] .+= Λ2[g,q].*view(Xc,q,:)
            RHS = view(Xc,q,:)'*view(yCorr,g,:) .+ α*pMeans[g]
            LHS = view(Xc,q,:)'*view(Xc,q,:)
            meanBeta = RHS/LHS
            nowBeta = sampleBeta(meanBeta, LHS, σ2ϵ)
            Λ2[g,q] = nowBeta
            yCorr[g,:] .-= nowBeta.*view(Xc,q,:)
        end
    end
end

#function forming Mj matrix.... Never changes
function lambda1BIGM(nG,Y,N)
    #sorted traits within individuals. nGenes lines per individual
    #Will be used for lambda1
    BIGY = Array{Array{Float64,2},1}()
    for i in 1:N
        local bigY = Y[2:end,i]'
    for g in 2:nG
        nowGenes = deleteat!(collect(1:nG), g)
        bigY = cat(bigY, Y[nowGenes,i]';dims=(1,2))
    end
        BIGY = push!(BIGY,bigY)
    end
    return BIGY #should it be sparse?
end

function logTargetDist(Λ,λ,λt,iVλ,N)
     return ((N/2)*log(det(Λ))) + (-0.5*((λ-λt)'*iVλ*(λ-λt)))
end
                        
function proposalDist(meanVec,varCov)
    return rand(MvNormal(meanVec,Symmetric(varCov)))
end

            
#main function to sample Lambda1
function sampleΛ1!(yCorr,N,MpM,BIGM,posΛ1,nλ1Coeff,MpM_I,starI,currentλ1,currentΛ1,accept,proposedΛ1;λ1_0=0.0,σ2ϵ=0.05,τ2=0.005)
    #add effect
    for i in eachindex(BIGM)
        yCorr[:,i] .+= BIGM[i]*currentΛ1[posΛ1]
    end
                                
    ### M-H here
    iLHS    = inv(MpM + MpM_I*(σ2ϵ/τ2))
    RHS     = ones(nλ1Coeff)*λ1_0*(σ2ϵ/τ2)
    for i in eachindex(BIGM)
        RHS .+= BIGM[i]'*yCorr[:,i]
    end
    
    λ1_mean = iLHS*RHS
    λ1_coVar  = iLHS*σ2ϵ

    #dont know what to do with these
    λt = λ1_mean        
    iλt_coVar = inv(λ1_coVar)
    ####
        
    proposedλ1 = proposalDist(λ1_mean,λ1_coVar)
    proposedΛ1[posΛ1] .= proposedλ1
    
    proposedΛstar = starI-proposedΛ1
    currentΛstar  = starI-currentΛ1
 
    A = exp(logTargetDist(proposedΛstar,proposedλ1,λt,iλt_coVar,N) - logTargetDist(currentΛstar,currentλ1,λt,iλt_coVar,N))
    if rand() < A
        currentλ1 .= proposedλ1       # accept move with probabily min(1,A)
        currentΛ1 .= proposedΛ1       # accept move with probabily min(1,A)
        accept  .+= 1                  # otherwise "reject" move, and stay where we are
    end

    ### 
                
    #remove effect
    for i in eachindex(BIGM)
        yCorr[:,i] .-= BIGM[i]*currentΛ1[posΛ1]
    end
end

function sampleBeta(meanBeta, lhs, σ2E)
    return rand(Normal(meanBeta,sqrt(lhs\σ2E)))
end

function sampleVarLambda1(νS,dataVec,df,Kstar)
    return((νS + dot(dataVec,dataVec))/rand(Chisq(df + Kstar)))
end
            
function sampleVarBeta(νS_β,whichLoci,df_β,regionSize)
    return((νS_β + dot(whichLoci,whichLoci))/rand(Chisq(df_β + regionSize)))
end
function sampleVarE(νS_e,yCorVec,df_e,nRecords)
    return((νS_e + dot(yCorVec,yCorVec))/rand(Chisq(df_e + nRecords)))
end
            
function outMCMC(folder::String,thisVar,output)
        out0 = open(folder*"/$(thisVar)Out", "a")
        writedlm(out0, output)
        close(out0)
end
            
function summaryMCMC(param;summary=false,outFolder=pwd()*"/outMCMC")
        param = CSV.read("$outFolder/$(param)Out",CSV.Tables.matrix,header=false)
                if summary==true
                        chn = Chains(param)
                        display(chn)
                        display(plot(chn))
                        param = mean(Matrix(param),dims=1)
                else param = mean(Matrix(param),dims=1)
                end
        return param
end
            
folderHandler = function(outFolder)
    if isdir(outFolder)==true
        println("Output folder $outFolder exists. Removing it")
        rm(outFolder,force=true,recursive=true)
        else
        println("$outFolder has been created to store the MCMC output")
        run(`mkdir $outFolder`)
    end
end


end
