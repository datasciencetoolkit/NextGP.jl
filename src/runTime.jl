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

struct BayesPRType
    r::Int
    m::Union{Vector{Float64},Float64}
    v::Union{Matrix{Float64},Float64}
    name::String
end

BayesPR(r::Int,m::Union{Vector{Float64},Float64},v::Union{Matrix{Float64},Float64};name="BayesPR") = BayesPRType(r,m,v,name)


struct RandomEffectType <: AbstractTerm
    str::Any
    m::Any
    v::Float64
end

Random(str::Any,m::Any,v::Float64) = RandomEffectType(str,m,v)

