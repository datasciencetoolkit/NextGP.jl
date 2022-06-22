module IO

export outMCMC
export summaryMCMC

using DelimitedFiles
using CSV
using StatsBase
using MCMCChains

macro name(arg)
    x = string(arg)
    quote
        $x
    end
end

function outMCMC(folder::String,thisVar,output)
        out0 = open(pwd()*"/$(thisVar)Out", "a")
        writedlm(out0, output)
        close(out0)
end

function summaryMCMC(param;summary=false)
	param = CSV.read(param*"Out",CSV.Tables.matrix,header=false)
		if summary==true
			chn = Chains(param)
			display(chn)
			plot(chn)
			param = mean(Matrix(param),dims=1)
		else param = mean(Matrix(param),dims=1)
		end
	return param	
end

end
