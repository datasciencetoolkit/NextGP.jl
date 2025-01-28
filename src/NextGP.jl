module NextGP

#exporting run-time equivalent of functions
export PED
export SNP
export BayesPR,BayesB,BayesC,BayesR,BayesRCπ,BayesRCplus,BayesLV
export Random
export BayesPRType
export SummaryStatistics


#
export makeA
export makePed
export makeG
export summaryMCMC
export show #import is in misc
export @model

using DataFrames
using CategoricalArrays
using Printf

include("misc.jl")
include("prepMatVec.jl")
include("types.jl")
include("MCMC.jl")
include("outFiles.jl")
include("GRN.jl")

using .MCMC
export runLMEM
using .prepMatVec
export prep
using .GRN
export estGRN_MHGibbs 


end
