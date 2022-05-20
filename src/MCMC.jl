module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers
using StatsBase

include("equations.jl")
include("runTime.jl")
include("samplers.jl")

struct MME
	Y
	X::Vector{Matrix{Float64}}
	Z::Vector{Matrix{Float64}}
        namesX
	namesZ	
end	

runGibbs = function(formula, userHints, userData, userPedData, chainLength, blockThese)
	yVec,FE,RE,namesFE,namesRE = equations.mme(formula, userHints, userData, userPedData, blockThese)
        varResidual = 35  ##### FIXED for now
        samplers.runSampler(yVec,FE,RE,varResidual,chainLength)
        return(yVec,FE,RE,namesFE)
end

end
