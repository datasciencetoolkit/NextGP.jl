module IO

export outMCMC

using DelimitedFiles

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

function getResults(param;summary=false)
	param = CSV.read(param*"Out",Tables.matrix,header=false)
		if summary=true
			chn = Chains(param)
			param = mean(Matrix(param),dims=1)
		else param = mean(Matrix(param),dims=1)
	return param	
end

end
