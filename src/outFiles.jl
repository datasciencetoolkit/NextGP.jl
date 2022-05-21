module IO

export outMCMC

using DelimitedFiles

function outMCMC(folder::String,b)
        out0 = open(pwd()*"/bOut", "a")
        writedlm(out0, b)
        close(out0)
end

end
