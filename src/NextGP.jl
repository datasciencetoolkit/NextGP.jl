module NextGP

#exporting run-time equivalent of functions
export ran

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels

include("MME.jl")
include("runTime.jl")
#include("MCMC.jl")

using .MCMC
export runGibbs

end
