#make design matrix for ran()

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

#
ran(arg1,arg2) = make_ran_matrix(arg2[!,Symbol(arg1)],arg2[!,Symbol(arg1)])
