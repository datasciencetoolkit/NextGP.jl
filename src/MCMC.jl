module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers
using StatsBase
using Printf

####REMOVE
using InteractiveUtils
####


include("equations.jl")
include("runTime.jl")
include("samplers.jl")


runGibbs = function(formula,userHints,userData,userPedData,nChain,nBurn,nThin,blockThese,VCV;genotypes...)
	idY,yVec,FE,RE,ME,regionSizes,namesFE,namesRE,namesME = equations.mme(formula,userHints,userData,userPedData,blockThese;paths2geno=genotypes)
        samplers.runSampler(idY,yVec,FE,RE,nChain,nBurn,nThin,VCV,ME,regionSizes)
	println("region sizes: $regionSizes")
        return(yVec,FE,RE,ME,namesFE,namesRE,namesME)
end

end
