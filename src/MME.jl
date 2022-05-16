module MME

export mme,ran

function make_ran_matrix(x1::AbstractVector,x2::AbstractVector)
           isa(x1|x2, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(x)))"))
           u = unique(x2)
           Z = Matrix{Int32}(undef, length(x1), length(u)) #Matrix{Bool}
           for i in eachindex(u)
               @. Z[:, i] = x1 .== u[i]
           end
           return Z
       end

ran(arg1,arg2) = make_ran_matrix(arg2[!,Symbol(arg1)])

function mme(f, userHints, userData)

	terms4StatsModels = String.(split(repr(f.rhs), ('+')))
	terms4StatsModels = replace.(terms4StatsModels, ":" => "")
	terms4StatsModels = [filter(x -> !isspace(x), trm) for trm in terms4StatsModels]

	yVec = StatsModels.modelmatrix(f.lhs, userData)

	FE = Array{Array{Float64,2},1}(undef,0)
	namesFE = []

	RE = Array{Array{Float64,2},1}(undef,0)
	namesRE = []

	for i in 1:length(f.rhs)
		if f.rhs[i] isa FunctionTerm{typeof(ran)}
			arg1 = repr((f.rhs[i].args_parsed)[1]) #now it is Symbol
			arg2 = repr((f.rhs[i].args_parsed)[2]) #now it is from string
                	arg2 = eval(Meta.parse(arg2)) #now it is from string to data. Later will be path

			println("$i has type ran Type")
                	println(ran(arg1, arg2))

                	push!(RE,ran(arg1, arg2))
                	push!(namesRE, terms4StatsModels[i])
		elseif f.rhs[i] isa FunctionTerm{typeof(|)} #to avoid schema issues/errors
			println("$i has type | Type")
			my_sch = schema(userData, userHints) #work on userData and userHints
			my_ApplySch = apply_schema(terms(f.rhs[i]), my_sch, MixedModel)
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
