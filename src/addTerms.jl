###make design matrices

using CategoricalArrays

function make_ran_matrix(x1::AbstractVector,x2::AbstractVector)
           isa(x1, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))
	   isa(x2, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))
           u = unique(x2);filter!(x->x≠0,u)
           Z = Matrix{Bool}(undef, length(x1), length(u))
           for i in eachindex(u)
               @. Z[:, i] = x1 .== u[i]
           end
           return Z
       end

ranMat(arg1,arg2,data1,data2) = make_ran_matrix(data1[!,Symbol(arg1)],data2[!,Symbol(arg2)])

