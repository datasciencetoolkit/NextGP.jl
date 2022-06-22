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


runGibbs = function(formula,userData,nChain,nBurn,nThin;userHints=Dict{String,Any}(),blockThese=[],VCV=[],userPedData=[],map=[],genotypes...)
	idY,A,yVec,FE,RE,ME,regionSizes = equations.mme(formula,userHints,userData,blocks=blockThese,path2ped=userPedData,paths2geno=genotypes)
        samplers.runSampler(idY,A,yVec,FE,RE,nChain,nBurn,nThin,VCV,ME,map,regionSizes)
        return(yVec,FE,RE,ME)
end

end
