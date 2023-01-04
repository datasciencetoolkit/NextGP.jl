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
include("mme.jl")
include("samplers.jl")
include("misc.jl")
include("outFiles.jl")

runLMEM = function(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],summaryStat=Dict{Any,Any}())
	
	folderHandler(outFolder)

	yVec,X,Z,M = prepMatVec.prep(formula,userData,userHints=myHints,path2ped=userPedData,priorVCV=VCV)

	ycorr,nData,dfE,scaleE,X,b,Z,u,varU,M,beta,varBeta = mme.getMME!(yVec,X,Z,M,blockThese,VCV,summaryStat,outFolder)
	
	samplers.runSampler!(ycorr,nData,dfE,scaleE,X,b,Z,u,varU,M,beta,varBeta,nChain,nBurn,nThin,outFolder)
	
end

end
