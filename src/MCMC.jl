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
include("misc.jl")
include("outFiles.jl")

runGibbs = function(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],map=[],genotypes...)
	
	folderHandler(outFolder)

	levelsFR,Ainv,yVec,FE,RE,ME,regionSizes = equations.mme(formula,userData,userHints=myHints,path2ped=userPedData,paths2geno=genotypes)
	
        samplers.runSampler(Ainv,yVec,FE,RE,levelsFR,blockThese,nChain,nBurn,nThin,VCV,ME,map,regionSizes,outFolder)

	#return(yVec,FE,RE,ME)
end

end
