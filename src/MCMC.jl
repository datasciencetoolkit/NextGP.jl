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


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese,VCV;genotypes...)
	idY,yVec,FE,RE,namesFE,namesRE = equations.mme(formula,userHints,userData,userPedData,blockThese;genotypes)
        samplers.runSampler(idY,yVec,FE,RE,nChain,nBurn,nThin,VCV)
        return(yVec,FE,RE,namesFE,namesRE)
end

end
