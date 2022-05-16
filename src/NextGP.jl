module NextGP

using DataFrames,CategoricalArrays,StatsModels,MixedModels
include("mme.jl")

runGibbs = function(formula::FormulaTerm, userHints::Dict, userData::DataFrame)
	return mme(formula, userHints, userData)
end

end
