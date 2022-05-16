module NextGP

#export runGibbs
export MME,MME.mme,MME.ran

import DataFrames
import CategoricalArrays
import StatsModels
import MixedModels

include("MME.jl")

runGibbs = function(formula, userHints, userData)
	return MME.mme(formula, userHints, userData)
end

end
