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


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese,iVCV,priorVCV)
	idY,yVec,FE,RE,namesFE,namesRE = equations.mme(formula,userHints,userData,userPedData,blockThese)
        samplers.runSampler(idY,yVec,FE,RE,nChain,nBurn,nThin,iVCV,priorVCV)
        return(yVec,FE,RE,namesFE,namesRE)
end

end
