###make design matrices

using CategoricalArrays

function make_ran_matrix(x1::AbstractVector,x2::AbstractVector)
           isa(x1, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))
	   isa(x2, CategoricalArray) ||
                       throw(ArgumentError("ran() only works with CategoricalArrays (got $(typeof(2)))"))
           
           filter!(x->x≠0,x1);filter!(x->x≠0,u2)
           if length(unique(x1)) < length(unique(x2))
		u2 = unique(x2)
                u1 = unique(x1)
           else
                u2 = unique(x1)
                u1 = unique(x2)
           end
	  # u = unique(x2);
          # filter!(x->x≠0,u)
          # Z = Matrix{Bool}(undef, length(x1), length(u))
           Z = Matrix{Bool}(undef, length(u1), length(u2))
           for i in eachindex(u)
               @. Z[:, i] = x1 .== u[i]
           end
           return Z
       end

ranMat(arg1,arg2,data1,data2) = make_ran_matrix(data1[!,Symbol(arg1)],data2[!,Symbol(arg2)])

