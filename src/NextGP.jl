module NextGP

#exporting run-time equivalent of functions
export ran
export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels

include("MME.jl")

runGibbs = function(formula, userHints, userData)
	return mme(formula, userHints, userData)
end

end
