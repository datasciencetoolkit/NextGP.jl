module NextGP

#exporting run-time equivalent of functions
export ran
export PR

#
export makeA
export makePed

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Printf

include("prepMatVec.jl")
include("runTime.jl")
include("MCMC.jl")
include("misc.jl")
include("outFiles.jl")

using .MCMC
export runLMEM
using .IO
export summaryMCMC
using .prepMatVec
export prep

end
