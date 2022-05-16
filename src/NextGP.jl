module NextGP

export ran,runGibbs

import DataFrames
import CategoricalArrays
import StatsModels
import MixedModels

include("MME.jl")
include("addTerms.jl")

runGibbs = function(formula, userHints, userData)
	return MME.mme(formula, userHints, userData)
end

end
