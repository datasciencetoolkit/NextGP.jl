

struct RandTerm <: AbstractTerm
    term::Symbol
    data::DataFrame
end

ran(s::Symbol, d::DataFrame) = RandTerm(term(s), d)

