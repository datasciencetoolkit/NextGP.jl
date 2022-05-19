module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers
using StatsBase

include("MME.jl")
include("runTime.jl")
include("samplers.jl")


runGibbs = function(formula, userHints, userData, userPedData, chainLength, blockThese)
	yVec,FE,RE,namesFE,namesRE = MME.mme(formula, userHints, userData, userPedData, blockThese)
        varResidual = 35  ##### FIXED for now
        samplers.runSampler(yVec,FE,RE,varResidual,chainLength)
        return(yVec,FE,RE,namesFE)
end

end
