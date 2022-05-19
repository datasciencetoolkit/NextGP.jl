module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers

include("MME.jl")
include("runTime.jl")
include("samplers.jl")


runGibbs = function(formula, userHints, userData, userPedData,chainLength)
	yVec,FE,RE,namesFE,namesRE = MME.mme(formula, userHints, userData, userPedData)
        varResidual = 350  ##### FIXED for now
        return  samplers.runSampler(yVec,FE,RE,varResidual,chainLength)
end

end
