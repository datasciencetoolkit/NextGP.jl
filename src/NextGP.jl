module NextGP

using DataFrames,CategoricalArrays,StatsModels,MixedModels
iclude("MME.jl")

runGibbs = function(formula::FormulaTerm, userHints::Dict, userData::DataFrame)
	return MME.mme(formula, userHints, userData)
end

end
