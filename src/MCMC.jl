module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels

include("MME.jl")
include("runTime.jl")

runGibbs = function(formula, userHints, userData, userPedData)
	return MME.mme(formula, userHints, userData, userPedData)
end

end
