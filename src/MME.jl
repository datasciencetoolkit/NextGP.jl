module MME

using StatsModels, MixedModels, CategoricalArrays

function make_ran_matrix(x1::AbstractVector,x2::AbstractVector)
           isa(x1, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))
           isa(x2, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))

           u = unique(x2);
           filter!(x->x≠0,u)
           Z = Matrix{Bool}(undef, length(x1), length(u))
           for i in eachindex(u)
               @. Z[:, i] = x1 .== u[i]
           end
           return Z
       end

ranMat(arg1,arg2,data1,data2) = make_ran_matrix(data1[!,Symbol(arg1)],data2[!,Symbol(arg2)])


function mme(f, userHints, userData, userPedData, blocks)
        terms4StatsModels = String.(split(repr(f.rhs), ('+')))
        terms4StatsModels = replace.(terms4StatsModels, ":" => "")
        terms4StatsModels = [filter(x -> !isspace(x), trm) for trm in terms4StatsModels]

        yVec = StatsModels.modelmatrix(f.lhs, userData)

        FE = Array{Array{Float64,2},1}(undef,0)
        namesFE = []

        RE = Array{Array{Float64,2},1}(undef,0)
        namesRE = []

        for i in 1:length(f.rhs)
                if (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "ran")
                        println("$i has type ran Type")
                        sym1 = repr((f.rhs[i].args_parsed)[1]) #now it is Symbol
                        sym2 = repr((f.rhs[i].args_parsed)[2]) #now it is from string
#                       arg2 = eval(Meta.parse(arg2)) #now it is from string to data. Later will be path
                        println("sym1: $sym1 sym2: $sym2")

                        push!(RE,ranMat(sym1, sym2, userData, userPedData))
                        push!(namesRE, terms4StatsModels[i])
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "|")
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

	mergedOnes = []
	delThese = []
	newNamesFE = []
	for i in blocks
        	println("blocking: $i")
		blockThese = findall(x->x in i, namesFE)
        	println("blocking pos: $blockThese")
		mergeTo = minimum(blockThese)
        	println("merged old pos: $mergeTo")
		FE[mergeTo] = hcat(FE[blockThese]...)
		delThese = vcat(delThese,blockThese[blockThese .!== minimum(blockThese)])
		push!(mergedOnes, blockThese)
	end

	deleteat!(outFE, sort(delThese))
        
        return yVec, FE, RE, namesFE, namesRE
        end

end

