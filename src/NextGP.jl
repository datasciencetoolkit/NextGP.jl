module NextGP

export BayesPR,BayesB,BayesC,BayesR,BayesRCπ,BayesRCplus,BayesLV

#
export makeA
export makePed
export makeG
export summaryMCMC
export show #import is in misc

using DataFrames
using CategoricalArrays
using Printf

include("misc.jl")

include("prepMatVec.jl")
#using .prepMatVec
#exporting run-time equivalent of functions
export @model
export Random,PED,SNP,BayesPRType,SummaryStatistics,
export DataTerm,ConstantTerm,FixedEffect

include("MCMC.jl")
include("outFiles.jl")
include("GRN.jl")

#using .MCMC
export runLMEM
using .GRN
export estGRN_MHGibbs 

end
