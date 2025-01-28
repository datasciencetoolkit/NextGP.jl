module MCMC

export runLMEM

using DataFrames
using CategoricalArrays
using Distributions,LinearAlgebra #samplers
using StatsBase
using Printf



include("prepMatVec.jl")
include("mme.jl")
include("samplers.jl")
#include("misc.jl")
include("outFiles.jl")


"""
	function runLMEM(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],summaryStat=Dict{Any,Any}())

* `formula` is the model equatio as made available by the [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package
* `userData` is a `DataFrame` including `lhs` and `rhs` variables (other than defined by `PED` and `SNP`)
* `nChain`,`nBurn` and `nThin` are the chain length, burn-in period, and the thining interval used for the MCMC sampling
* Users can define coding of their variables (e.g. full dummy coding) by providing `myHints`. Check `StatsModels.jl`'s manual for [`categorical data`](https://juliastats.org/StatsModels.jl/latest/contrasts/#Modeling-categorical-data) could be useful.  
"""
runLMEM = function(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],summaryStat=Dict{Any,Any}())
	
	folderHandler(outFolder)

	yVec,X,Z,M = prepMatVec.prep(formula,userData,userHints=myHints,path2ped=userPedData,priorVCV=VCV)

	ycorr,nData,E,X,b,Z,u,varU,M,beta,varBeta,delta = mme.getMME!(yVec,X,Z,M,blockThese,VCV,summaryStat,outFolder)

	samplers.runSampler!(ycorr,nData,E,X,b,Z,u,varU,M,beta,varBeta,delta,nChain,nBurn,nThin,outFolder)
	
end

end
