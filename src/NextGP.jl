module NextGP

using DataFrames,CategoricalArrays,StatsModels,MixedModels
include("MME.jl")

runGibbs = function(formula, userHints, userData)
	return MME.mme(formula, userHints, userData)
end

end
