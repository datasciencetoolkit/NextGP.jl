module NextGP

#exporting run-time equivalent of functions
export ran
export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels

include("MME.jl")
include("runTime.jl")

import .RUNTIME: ran

export ran

runGibbs = function(formula, userHints, userData)
	return MME.mme(formula, userHints, userData)
end

end
