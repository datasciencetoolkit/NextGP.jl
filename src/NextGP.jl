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
using Distributions 
using LinearAlgebra
using StatsBase
using Printf
using CSV
using DataStructures
using PrettyTables
using ProgressMeter

include("types.jl")
#exporting run-time equivalent of functions
#export Random,PED,SNP,BayesPRType,SummaryStatistics
export Random,PED,SNP,GBLUP,SummaryStatistics

include("misc.jl")
include("outFiles.jl")
include("functions.jl")
using .functions
include("model.jl")
export @model
include("prepMatVec.jl")
include("mme.jl")
include("samplers.jl")
include("MCMC.jl")

include("GRN.jl")

export runLMEM

using .GRN
export estGRN_MHGibbs 

end
