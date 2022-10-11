using StatsModels: AbstractTerm

struct RandTerm <: AbstractTerm
    var::Char
    data::Char
end

ran(s::Char) = RandTerm(s)

struct PartRegTerm <: AbstractTerm
    mat::Char
    regSize::Int
    path::String
end

PR(s::Char, d::Int, p::String) = PartRegTerm(s, d, p)

