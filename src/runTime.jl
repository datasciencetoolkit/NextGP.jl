module RUNTIME

using StatsModels: AbstractTerm
using DataFrames

struct RandTerm <: AbstractTerm
    term::Symbol
    data::DataFrame
end

ran(s::Symbol, d::DataFrame) = RandTerm(term(s), d)

end
