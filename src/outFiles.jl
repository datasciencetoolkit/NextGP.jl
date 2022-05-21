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
	varName = @name thisVar
        out0 = open(pwd()*"/$(thisVar)Out", "a")
        writedlm(out0, output)
        close(out0)
end

end
