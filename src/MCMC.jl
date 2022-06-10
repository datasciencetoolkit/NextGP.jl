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


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese,VCV,VCVM;map,genotypes...)
	idY,yVec,FE,RE,ME,regionSizes = equations.mme(formula,userHints,userData,userPedData,blockThese;paths2geno=genotypes)
        samplers.runSampler(idY,yVec,FE,RE,nChain,nBurn,nThin,VCV,VCVM,ME,map,regionSizes)
	println("region sizes: $regionSizes")
        return(yVec,FE,RE,ME)
end

end
