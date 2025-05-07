#module MCMC

#export runLMEM

#from prepMatVec
using CategoricalArrays, CSV, StatsBase, DataStructures, DataFrames, PrettyTables, LinearAlgebra

include("mme.jl")
include("samplers.jl")


"""
	function runLMEM(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],summaryStat=Dict{Any,Any}())

* `formula` is the model equatio as made available by the [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package
* `nChain`,`nBurn` and `nThin` are the chain length, burn-in period, and the thining interval used for the MCMC sampling
* Users can define coding of their variables (e.g. full dummy coding) by providing `myHints`. Check `StatsModels.jl`'s manual for [`categorical data`](https://juliastats.org/StatsModels.jl/latest/contrasts/#Modeling-categorical-data) could be useful.  
"""
runLMEM = function(model...;nChain=10000,nBurn=1000,nThin=10,myHints=Dict{Symbol,Any}(),blockThese=Dict{Symbol,Any}(),outFolder="outMCMC",VCV=Dict{Union{Expr,Symbol,Tuple}, Any}(),userPedData=[],summaryStat=Dict{Any,Any}())

	if length(model)==1 && isa(model[1].lhs,Symbol)
		println("I AM A lmm")
	elseif length(model)==1 && isa(model[1].lhs,Expr)
		println("I AM A MV lmm")
	elseif length(model)>1
		println("I AM A multi population lmm")
	else nothing
	end

	isa(modelType(model),Tuple{Vararg{lmm}}) ? println("I AM A lmm TYPE") : println("I AM NOT A lmm TYPE")

	folderHandler(outFolder)

	#Y,X,Z,M,E,modelInformation = prepMatVec.prep(model,path2ped=userPedData,priorVCV=VCV)
	Y,X,Z,M,E,modelInformation = prep(model,path2ped=userPedData,priorVCV=VCV)

	modelInformation,ycorr,nData,E,varE,X,b,Z,u,varU,M,beta,varBeta,delta = mme.getMME!(Y,X,Z,M,E,blockThese,VCV,summaryStat,modelInformation,outFolder)

	runSampler!(modelInformation,ycorr,nData,E,varE,X,b,Z,u,varU,M,beta,varBeta,delta,nChain,nBurn,nThin,outFolder)
	
end

#end
