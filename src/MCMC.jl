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

	println("Building parts of MME")

	levelsRE,levelsFE,Ainv,yVec,FE,RE,ME,regionSizes = equations.mme(formula,userData,userHints=myHints,blocks=blockThese,path2ped=userPedData,paths2geno=genotypes)


	#make MCMC output files
	IO.outMCMC(outFolder,"b",levelsFE)
	IO.outMCMC(outFolder,"u",levelsRE)

	IO.outMCMC(outFolder,"varE","varE")

	###
	levelsY = levelsRE	
	
        samplers.runSampler(levelsY,Ainv,yVec,FE,RE,nChain,nBurn,nThin,VCV,ME,map,regionSizes,outFolder)

	#return(yVec,FE,RE,ME)
end

end
