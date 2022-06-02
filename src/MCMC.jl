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


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese,VCV;map,genotypes...)
	idY,yVec,FE,RE,ME,regionSizes,namesFE,namesRE,namesME = equations.mme(formula,userHints,userData,userPedData,blockThese;paths2geno=genotypes)
        samplers.runSampler(idY,yVec,FE,RE,nChain,nBurn,nThin,VCV,ME,map,regionSizes)
	println("region sizes: $regionSizes")
        return(yVec,FE,RE,ME,namesFE,namesRE,namesME)
end

end
