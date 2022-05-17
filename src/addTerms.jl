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




#ran(arg1,arg2) = make_ran_matrix(arg2[!,Symbol(arg1)])

struct RandTerm <: AbstractTerm
    term::Symbol
    data::DataFrame
end

ran(s::Symbol, d::DataFrame) = RandTerm(term(s), d)

#ran = make_ran_matrix(ranMat.data[!,Symbol(ranMat.term)])
