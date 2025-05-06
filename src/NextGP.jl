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
export lmm

using DataFrames
using CategoricalArrays
using Printf

include("misc.jl")
include("prepMatVec.jl")
include("MCMC.jl")
include("outFiles.jl")
include("GRN.jl")
include("types.jl")

using .prepMatVec
export @model
export prep
using .MCMC
export runLMEM
using .GRN
export estGRN_MHGibbs 

export DataTerm


end
