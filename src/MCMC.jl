module MCMC

export runGibbs

using DataFrames
using CategoricalArrays
using StatsModels
using MixedModels
using Distributions,LinearAlgebra #samplers
using StatsBase
using Printf



include("equations.jl")
include("runTime.jl")
include("samplers.jl")
include("misc.jl")
include("outFiles.jl")

runGibbs = function(formula,userData,nChain,nBurn,nThin;myHints=Dict{Symbol,Any}(),blockThese=[],outFolder="outMCMC",VCV=[],userPedData=[],map=[],genotypes...)
	
	folderHandler(outFolder)

	println("Building parts of MME")

	levelsRE,levelsFE,Ainv,yVec,FE,RE,ME,regionSizes = equations.mme(formula,userData,userHints=myHints,blocks=blockThese,path2ped=userPedData,paths2geno=genotypes)


	#########make MCMC output files. Maybe better to move into samplers
	IO.outMCMC(outFolder,"b",levelsFE)

        for i in 1:length(levelsRE)
		nameRE = hcat(vcat(collect(values(levelsRE))[i]...)...)
		println("names for $(collect(keys(levelsRE))[i]): $(collect(values(levelsRE))[i])")
		IO.outMCMC(outFolder,"u$i",nameRE)
		IO.outMCMC(outFolder,"varU$i",join(collect(keys(levelsRE))[i],"_"))
	end	
	

	IO.outMCMC(outFolder,"varE",["varE"])
	##########
	
        samplers.runSampler(Ainv,yVec,FE,RE,nChain,nBurn,nThin,VCV,ME,map,regionSizes,outFolder)

	#return(yVec,FE,RE,ME)
end

end
