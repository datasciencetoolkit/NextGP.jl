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


runGibbs = function(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],VCV=[],userPedData=[],map=[],genotypes...)
	println("Building parts of MME")
	idY,A,yVec,FE,RE,ME,regionSizes = equations.mme(formula,userData,userHints=myHints,blocks=blockThese,path2ped=userPedData,paths2geno=genotypes)
	println("Running MCMC")
        samplers.runSampler(idY,A,yVec,FE,RE,nChain,nBurn,nThin,VCV,ME,map,regionSizes)
	#return(yVec,FE,RE,ME)
end

end
