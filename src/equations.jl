module equations

using StatsModels, MixedModels, CategoricalArrays, CSV, StatsBase, DataStructures

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
           return unique(x1),Z
       end

ranMat(arg1,arg2,data1,data2) = make_ran_matrix(data1[!,Symbol(arg1)],data2[!,Symbol(arg2)])


function mme(f, userHints, userData, userPedData, blocks; paths2geno)
        terms4StatsModels = String.(split(repr(f.rhs), ('+')))
        terms4StatsModels = replace.(terms4StatsModels, ":" => "")
        terms4StatsModels = [filter(x -> !isspace(x), trm) for trm in terms4StatsModels]

        yVec = StatsModels.modelmatrix(f.lhs, userData)
	
        FE = OrderedDict{Any,Any}() #any to block work

        RE = OrderedDict{String,Any}()

	ME = OrderedDict{String,Array{Float64, 2}}()
	regionSizes = Dict{String,Int64}()
		
	#column id within pedigree
	idRE = []
	
        for i in 1:length(f.rhs)
		if (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "PR")
			println("$i has type BayesPR Type")
			println("terms4StatsModels[i]: $(terms4StatsModels[i])")
			arg1 = repr((f.rhs[i].args_parsed)[1])
			arg2 = parse(Int64,repr((f.rhs[i].args_parsed)[2]))
			path = paths2geno[Symbol(arg1)]
			thisM = CSV.read(path,CSV.Tables.matrix)
			#centering
			thisM .-= mean(thisM,dims=1) 
			println("size of $arg1 data: $(size(thisM))")
			println("region size for $arg1: $arg2")
			ME[arg1] = thisM
                        thisM = 0 #I can directly merge to dict above
			regionSizes[arg1] = arg2
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "ran")
                        println("$i has type ran Type")
                        sym1 = repr((f.rhs[i].args_parsed)[1]) #now it is Symbol
                        sym2 = repr((f.rhs[i].args_parsed)[2]) #now it is from string
#                       arg2 = eval(Meta.parse(arg2)) #now it is from string to data. Later will be path
                        println("sym1: $sym1 sym2: $sym2")
			
			IDs,thisZ = ranMat(sym1, sym2, userData, userPedData)
			RE[sym1] = thisZ
			thisZ = 0
			push!(idRE,IDs)
                elseif (f.rhs[i] isa FunctionTerm) && (String(nameof(f.rhs[i].forig)) == "|")
                        println("$i has type | Type")
                        my_sch = schema(userData, userHints) #work on userData and userHints
                        my_ApplySch = apply_schema(terms(f.rhs[i]), my_sch, MixedModels.MixedModel)
			#####NO IDs for this effect!!! Will be added later!!!!#####################################################
                       	thisZ = modelcols(my_ApplySch, userData)
			RE[terms4StatsModels[i]] = thisZ
			thisZ = 0
			push!(RE,modelcols(my_ApplySch, userData))

                else
                println("$i has type $(typeof(f.rhs[i]))")
		thisX = StatsModels.modelmatrix(f.rhs[i], userData,hints= userHints)
		FE[terms4StatsModels[i]] = thisX
		thisX = 0
                end
        end

	#BLOCK FIXED EFFECTS
	for b in blocks
	getThese = intersect(collect(keys(FE)), b)
	FE[Tuple(getThese)] = hcat(getindex.(Ref(FE), getThese)...)
	for d in getThese
		delete!(FE,d)
	end
end


        
        return idRE, vec(yVec), FE, RE, ME, regionSizes
        end

end

