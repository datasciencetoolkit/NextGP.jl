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


runGibbs = function(formula, userHints, userData, userPedData,chainLength)
	yVec,FE,RE,namesFE,namesRE = MME.mme(formula, userHints, userData, userPedData)
        varResidual = 350  ##### FIXED for now
        samplers.runSampler(yVec,FE,RE,varResidual,chainLength)
        return(yVec,FE,RE)
end

end
