module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers
using StatsBase
using Printf

include("equations.jl")
include("runTime.jl")
include("samplers.jl")


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese,iVCV)
	idY,yVec,FE,RE,namesFE,namesRE = equations.mme(formula,userHints,userData,userPedData,blockThese)
        varResidual = 350  ##### FIXED for now
	varRandomEffects = 150 ### 150 for a. FIXED for now
        samplers.runSampler(idY,yVec,FE,RE,varResidual,varRandomEffects,nChain,nBurn,nThin,iVCV)
        return(yVec,FE,RE,namesFE,namesRE)
end

end
