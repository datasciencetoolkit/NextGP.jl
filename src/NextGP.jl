module NextGP

#exporting run-time equivalent of functions
export PED
export SNP
export BayesPR
export Random
export BayesPRType
export SummaryStatistics
export estGRN_MHGibbs 


#
export makeA
export makePed
export makeG

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
