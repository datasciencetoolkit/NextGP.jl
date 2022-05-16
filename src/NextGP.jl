module NextGP

#exporting ran from MME
export ran,runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels

include("MME.jl")
include("addTerms.jl")

runGibbs = function(formula, userHints, userData)
	return MME.mme(formula, userHints, userData)
end

end
