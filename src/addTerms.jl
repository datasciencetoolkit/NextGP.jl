#make design matrix for ran()

using CategoricalArrays

function make_ran_matrix(x1::AbstractVector)
           isa(x1, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(x)))"))
           u = unique(x1)
           Z = Matrix{Bool}(undef, length(x1), length(u))
           for i in eachindex(u)
               @. Z[:, i] = x1 .== u[i]
           end
           return Z
       end




