module MCMC

export runLMEM

using DataFrames
using CategoricalArrays
using Distributions,LinearAlgebra #samplers
using StatsBase
using Printf

include("types.jl")
include("prepMatVec.jl")
include("mme.jl")
include("samplers.jl")
include("misc.jl")
include("outFiles.jl")


"""
	function runLMEM(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],summaryStat=Dict{Any,Any}())

* `formula` is the model equatio as made available by the [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package
* `nChain`,`nBurn` and `nThin` are the chain length, burn-in period, and the thining interval used for the MCMC sampling
* Users can define coding of their variables (e.g. full dummy coding) by providing `myHints`. Check `StatsModels.jl`'s manual for [`categorical data`](https://juliastats.org/StatsModels.jl/latest/contrasts/#Modeling-categorical-data) could be useful.  
"""
runLMEM = function(model...;nChain=10000,nBurn=1000,nThin=10,myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],summaryStat=Dict{Any,Any}())

	println("model is a $(typeof(model))")

	isa(typeof{model},Tuple{Vararg{lmm}}) || isa(typeof{model},Tuple{lmm}) ? println("ALL IS WELL") : throw(ArgumentError("Please enter a valid model"))
	
	isa(model,Tuple{Vararg{lmm}}) || isa(model,Tuple{lmm}) ? nothing : throw(ArgumentError("Please enter a valid model"))

	folderHandler(outFolder)

	yVec,X,Z,M = prepMatVec.prep(model,path2ped=userPedData,priorVCV=VCV)

	ycorr,nData,E,X,b,Z,u,varU,M,beta,varBeta,delta = mme.getMME!(yVec,X,Z,M,blockThese,VCV,summaryStat,outFolder)

	samplers.runSampler!(ycorr,nData,E,X,b,Z,u,varU,M,beta,varBeta,delta,nChain,nBurn,nThin,outFolder)
	
end

end
