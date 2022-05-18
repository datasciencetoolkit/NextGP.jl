using StatsModels: AbstractTerm

struct RandTerm <: AbstractTerm
    var::Char
    data::Char
end

ran(s::Char, d::Char) = RandTerm(s, d)

