module inOut

export outMCMC

using DelimitedFiles
using CSV
using StatsBase
using DataFrames

macro name(arg)
    x = string(arg)
    quote
        $x
    end
end

function outMCMC(folder::String,thisVar,output)
        out0 = open(folder*"/$(thisVar)Out", "a")
        writedlm(out0, output)
        close(out0)
end


end
