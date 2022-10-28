module MCMC

export runLMEM

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers
using StatsBase
using Printf



include("prepMatVec.jl")
include("runTime.jl")
include("samplers.jl")
include("misc.jl")
include("outFiles.jl")

runLMEM = function(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],map=[],genotypes...)
	
	folderHandler(outFolder)

	levelsFR,Ainv,yVec,FE,RE,ME,regionSizes = prepMatVec.prep(formula,userData,userHints=myHints,path2ped=userPedData,paths2geno=genotypes)
	
        samplers.runSampler(Ainv,yVec,FE,RE,levelsFR,blockThese,nChain,nBurn,nThin,VCV,ME,map,regionSizes,outFolder)

	#return(yVec,FE,RE,ME)
end

end
