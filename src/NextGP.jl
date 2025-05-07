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

#include("types.jl")
#exporting run-time equivalent of functions
#export Random,PED,SNP,BayesPRType,SummaryStatistics

include("misc.jl")
include("prepMatVec.jl")
include("MCMC.jl")
include("outFiles.jl")
include("GRN.jl")

include("model.jl")
export @model

using .prepMatVec
export Random,PED,SNP,BayesPRType,SummaryStatistics
export DataTerm


#export @model,prep

using .MCMC
export runLMEM
using .GRN
export estGRN_MHGibbs 

end
