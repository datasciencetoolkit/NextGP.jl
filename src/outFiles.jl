module IO

using DelimitedFiles

outMCMC = function(folder::String,b::Array{Float64})
        out0 = open(pwd()*"/muOut", "a")
        writedlm(out0, b)
        close(out0)
end

end
