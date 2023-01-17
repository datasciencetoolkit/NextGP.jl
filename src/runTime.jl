using StatsModels

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

"""
        SNP(name::Char,path::String;map::String="")
* Defines SNP information for further analysis.
* `path` is the path for the marker file.
* Marker files are currently expected to be ordered as the phenotype data.
* The method to be applied to the data, `GBLUP`, `BayesB`, `BayesR` etc, is defined in the prior setting.
* Map file is optional. If not provided, a Bayesian Regression model with common variance for all SNPs will be applied. If provided, shoul match the order in the genotype file.
* One most avoid overlapping marker sets by using different `name`s.
"""
SNP(name::Char,path::String;map::String="") = GenomicTerm(name,path,map)

struct BayesPRType
    r::Int
    v::Union{Matrix{Float64},Float64}
    name::String
end

"""
        BayesPR(r::Int,v::Union{Matrix{Float64},Float64};name="BayesPR")
* `r` is the region size. In other words, the number of SNPs that share a common variance.
    * `1`: each SNP has its own (co)variance
    * `99`: SNPs on the same chromosome has the same (co)variance
    * `9999`: All SNPs have the same (co)variance
    * One can define any other region size, for example, 30, 40 or 100.
* `v` is the variance for the prior distribution of SNPs.
"""
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
    class::Vector{Float64}
    v::VarCovarTypes
    name::String
end

BayesR(pi::PiTypes,class::Vector{Float64},v::VarCovarTypes;name="BayesR") = BayesRType(pi,class,v,name)

struct BayesLogVarType
    v::Union{Matrix{Float64},Float64}
    f::StatsModels.TermOrTerms
    covariates::DataFrame
    name::String
end

BayesLV(v::Float64,f::StatsModels.TermOrTerms,covariates::DataFrame;name="BayesLV") = BayesLogVarType(v,f,covariates,name)

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

