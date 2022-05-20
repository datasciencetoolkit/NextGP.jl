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


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese)
	yVec,FE,RE,namesFE,namesRE = equations.mme(formula,userHints,userData,userPedData,blockThese)
        varResidual = 350  ##### FIXED for now
        samplers.runSampler(yVec,FE,RE,varResidual,nChain,nBurn,nThin)
        return(yVec,FE,RE,namesFE,namesRE)
end

end
