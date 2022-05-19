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


runGibbs = function(formula, userHints, userData, userPedData)
	yVec,FE,RE,namesFE,namesRE = MME.mme(formula, userHints, userData, userPedData)
        varE ##### FIXED for now
        return  samplers.runSampler(Y,X,Z,varE)
end

