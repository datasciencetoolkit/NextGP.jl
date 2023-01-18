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
        function SNP(name::Char,path::String;map::String="")
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
        function BayesPR(r::Int,v::Union{Matrix{Float64},Float64})
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

"""
        function BayesB(pi::PiTypes,v::VarCovarTypes)
* `pi` is the proportion of SNPs to be included in the model at each MCMC cycel. 
* `v` is the variance for the prior distribution of SNPs.
"""
BayesB(pi::PiTypes,v::VarCovarTypes;name="BayesB") = BayesBType(pi,v,name)

struct BayesCType
    pi::PiTypes
    v::VarCovarTypes
    name::String
end

"""
        function BayesC(pi::PiTypes,v::VarCovarTypes)
* `pi` is the proportion of SNPs to be included in the model at each MCMC cycel. 
* `v` is the variance for the prior distribution of SNPs.
"""
BayesC(pi::PiTypes,v::VarCovarTypes;name="BayesC") = BayesCType(pi,v,name)

struct BayesRType
    pi::PiTypes
    class::Vector{Float64}
    v::VarCovarTypes
    name::String
end

"""
        function BayesR(pi::PiTypes,class::Vector{Float64},v::VarCovarTypes)
* `pi` is the vector of proportion of SNPs for each variance class.
* `class` is the vector of scales of common SPN variance for each variance class. The scales should be in the increasing order. For example, [0.0,0.0001,0.001,0.01].
* `v` is the variance for the prior distribution of SNPs.
"""
BayesR(pi::PiTypes,class::Vector{Float64},v::VarCovarTypes;name="BayesR") = BayesRType(pi,class,v,name)

struct BayesLogVarType
    v::Union{Matrix{Float64},Float64}
    f::StatsModels.TermOrTerms
    covariates::DataFrame
    name::String
end

"""
        function BayesLV(v::Float64,f::StatsModels.TermOrTerms,covariates::DataFrame)
* `v` is the variance for the prior distribution of SNPs.
* `f` is the model formula for the variance
* `covariates`is the `DataFrame` that includes explanatory varibles for the variance of each SNP. 
"""
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

