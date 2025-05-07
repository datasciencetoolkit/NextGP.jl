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
#exporting run-time equivalent of functions
export Random,PED,SNP,BayesPRType,SummaryStatistics,DataTerm,lmm,FixedEffect
export @model

include("MCMC.jl")
include("outFiles.jl")
include("GRN.jl")

using .MCMC
export runLMEM
using .GRN
export estGRN_MHGibbs 

end
