module MME

export FE,RE,namesFE,namesRE

include("addTerms.jl")

FE,RE,namesFE,namesRE = mme(f, userHints, userData, userPedData)

end
