#make design matrix for ran()

function make_ran_matrix(x::AbstractVector)
           isa(x, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(x)))"))
           u = unique(x)
           m = Matrix{Bool}(undef, length(x), length(u))
           for i in eachindex(u)
               @. m[:, i] = x .== u[i]
           end
           return m
       end


ran(arg1,arg2) = make_ran_matrix(arg2[!,Symbol(arg1)])

