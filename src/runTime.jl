
PiTypes = Union{Vector{Float64},Float64} #pi can be different (vector) or same per SNP (NO COR BayesB YET). BayesR also takes a vector of pi
VarCovarTypes = Union{Vector{Matrix{Float64}},Vector{Float64},Matrix{Float64},Float64} #prior for (co)var can be different (vector) or same per SNP


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
    pi::PiTypes
    v::VarCovarTypes
    name::String
end

BayesB(pi::PiTypes,v::VarCovarTypes;name="BayesB") = BayesBType(pi,v,name)

struct BayesCType
    pi::PiTypes
    v::VarCovarTypes
    name::String
end

BayesC(pi::PiTypes,v::VarCovarTypes;name="BayesC") = BayesCType(pi,v,name)

struct BayesRType
    pi::PiTypes
    v::VarCovarTypes
    name::String
end

BayesR(pi::PiTypes,v::VarCovarTypes;name="BayesR") = BayesRType(pi,v,name)


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

