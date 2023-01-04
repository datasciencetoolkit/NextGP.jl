
DiffOrSameπ = Union{Vector{Float64},Float64}
DiffOrSamePriorCoVar = Union{Vector{Matrix{Float64}},Vector{Float64},Matrix{Float64}


struct RandTerm
    var::Char
end

PED(s::Char) = RandTerm(s)

struct GenomicTerm
    name::Char
    path::String
    map::String
end

SNP(name::Char,path::String;map::String="") = GenomicTerm(name,path,map)

struct BayesPRType
    r::Int
    v::Union{Matrix{Float64},Float64}
    name::String
end

BayesPR(r::Int,v::Union{Matrix{Float64},Float64};name="BayesPR") = BayesPRType(r,v,name)


struct BayesBType
    π::DiffOrSameπ          #pi can be different (vector) or same per SNP
    v::DiffOrSamePriorCoVar #(co)var can be different (vector) or same per SNP
    name::String
end

BayesB(π::DiffOrSameπ,v::DiffOrSamePriorCoVar;name="BayesB") = BayesBType(π,v,name)


struct RandomEffectType
    str::Any
    v::Union{Matrix{Float64},Float64}
    type::Int
end

Random(str::Any,v::Union{Matrix{Float64},Float64};type=1) = RandomEffectType(str,v,type)


struct SummaryStatistics
    m::Union{Vector{Float64},Float64}
    v::Union{Vector{Float64},Matrix{Float64},Float64}
end

