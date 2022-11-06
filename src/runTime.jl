using StatsModels: AbstractTerm

struct RandTerm <: AbstractTerm
    var::Char
end

PED(s::Char) = RandTerm(s)

struct GenomicTerm <: AbstractTerm
    name::Char
    path::String
    map::String
end

SNP(name::Char,path::String;map::String="") = GenomicTerm(name,path,map)

struct BayesPRType <: AbstractTerm
    r::Int
    m::Any
    v::Float64
end

BayesPR(r::Int,m::Any,v::Float64) = BayesPRType(r,m,v)

struct RandomEffectType <: AbstractTerm
    str::Any
    m::Any
    v::Float64
end

Random(str::Any,m::Any,v::Float64) = RandomEffectType(str,m,v)

