module NextGP

#exporting run-time equivalent of functions
export ran
export PR

#
export makeA

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Printf

include("equations.jl")
include("runTime.jl")
include("MCMC.jl")
include("misc.jl")
include("outFiles.jl")

using .MCMC
export runGibbs
export summaryMCMC

end
