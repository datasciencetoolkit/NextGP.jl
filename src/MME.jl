module MME

include("addTerms.jl")

using StatsModels,CategoricalArrays,MixedModels

export mme

function mme(f, userHints, userData)
	terms4StatsModels = String.(split(repr(f.rhs), ('+')))
	terms4StatsModels = replace.(terms4StatsModels, ":" => "")
	terms4StatsModels = [filter(x -> !isspace(x), trm) for trm in terms4StatsModels]

	yVec = StatsModels.modelmatrix(f.lhs, userData)

	FE = Array{Array{Float64,2},1}(undef,0)
	namesFE = []

	RE = Array{Array{Float64,2},1}(undef,0)
	namesRE = []

println("data: $userData")
println("TYPE: $(f.rhs[5] isa FunctionTerm{typeof(ran)})")
println("$(typeof(f.rhs[5]))")
println("$(FunctionTerm{typeof(ran)})")

FunctionTerm{typeof(NextGP.MME.ran)} = FunctionTerm{typeof(ran)}

	for i in 1:length(f.rhs)
		if f.rhs[i] isa FunctionTerm{typeof(ran)}
			println("$i has type ran Type")			
			arg1 = repr((f.rhs[i].args_parsed)[1]) #now it is Symbol
			arg2 = repr((f.rhs[i].args_parsed)[2]) #now it is from string
                	arg2 = eval(Meta.parse(arg2)) #now it is from string to data. Later will be path
                        arg2 = userData 
			println("arg1: $arg1 arg2: $arg2")	
#                	println(ran(arg1, arg2))
		
#                	push!(RE,ran(arg1, arg2))
#                	push!(namesRE, terms4StatsModels[i])
		elseif f.rhs[i] isa FunctionTerm{typeof(|)} #to avoid schema issues/errors
			println("$i has type | Type")
			my_sch = schema(userData, userHints) #work on userData and userHints
			my_ApplySch = apply_schema(terms(f.rhs[i]), my_sch, MixedModels.MixedModel)
			println(modelcols(my_ApplySch, userData)) #work on userData and userHints
               	 	push!(RE,modelcols(my_ApplySch, userData))
               	 	push!(namesRE, terms4StatsModels[i])

		else
		println("$i has type $(typeof(f.rhs[i]))")
		println(StatsModels.modelmatrix(f.rhs[i], userData,hints= userHints)) #work on userData and userHints
                push!(FE,StatsModels.modelmatrix(f.rhs[i], userData,hints= userHints))
                push!(namesFE, terms4StatsModels[i])
		end
	end
	return FE, RE, namesFE, namesRE
	end
end
