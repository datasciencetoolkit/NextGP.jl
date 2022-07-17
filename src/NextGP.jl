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

include("equations.jl")
include("runTime.jl")
include("MCMC.jl")
include("misc.jl")
include("outFiles.jl")

using .MCMC
export runGibbs
using .IO
export summaryMCMC
using .equations

export mme

"""
    MyFunction(x,y)
This is an example of Docstring. This function receives two 
numbers x and y and returns the sum of the squares.
"""
function MyFunction(x,y)
    return x^2+y^2
end

export MyFunction

end
