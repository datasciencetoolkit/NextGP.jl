module NextGP

#exporting run-time equivalent of functions
export PED
export SNP
export BayesPR,BayesB,BayesC,BayesR
export Random
export BayesPRType
export SummaryStatistics


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
include("GRN.jl")

using .MCMC
export runLMEM
using .IO
export summaryMCMC
using .prepMatVec
export prep
using .GRN
export estGRN_MHGibbs 


end
