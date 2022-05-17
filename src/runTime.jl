using StatsModels: AbstractTerm
using DataFrames

struct RandTerm <: AbstractTerm
    var::Char
    data::DataFrame
end

ran(s::Char, d::DataFrame) = RandTerm(s, d)

