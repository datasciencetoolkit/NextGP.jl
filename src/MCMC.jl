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
	idY,yVec,FE,RE,GE,namesFE,namesRE,namesGE = equations.mme(formula,userHints,userData,userPedData,blockThese;paths2geno=genotypes)
	InteractiveUtils.varinfo()
        samplers.runSampler(idY,yVec,FE,RE,nChain,nBurn,nThin,VCV)
	InteractiveUtils.varinfo()
        return(yVec,FE,RE,GE,namesFE,namesRE,namesGE)
end

end
