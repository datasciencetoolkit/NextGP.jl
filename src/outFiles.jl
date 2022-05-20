module IO

export outMCMC

using DelimitedFiles

outMCMC = function(folder::String,b)
        out0 = open(pwd()*"/bOut", "a")
        writedlm(out0, b)
        close(out0)
end

end
