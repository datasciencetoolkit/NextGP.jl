using NextGP
using Test

@testset "NextGP.jl" begin
    	
	using StableRNGs; rng = StableRNG(1);
	data = DataFrame(y = rand(rng, 4), a = rand(rng, 4), b = [1:4;])
	data.b = CategoricalArray(data.b)

###################MODEL
f = @formula(y ~ 1 + a + poly(a,2) + (1|b) + ran(b,data)) #data will be path, and not phenotypic data

userHints = Dict(:b => StatsModels.FullDummyCoding())

print(runGibbs(formula,userHints,userData))

end
