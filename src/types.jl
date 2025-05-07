abstract type FixedEffect end


PiTypes = Union{Vector{Float64},Float64} #pi can be different (vector) or same per SNP (NO COR BayesB YET). BayesR also takes a vector of pi
VarCovarTypes = Union{Vector{Matrix{Float64}},Vector{Float64},Matrix{Float64},Float64} #prior for (co)var can be different (vector) or same per SNP
ExprOrSymbol = Union{Expr,Symbol}
ExprOrSymbolOrTuple = Union{Expr,Symbol,Tuple}

#the core lmm structure
struct lmm
	data::String
	model::Expr
	sym::Any
	lhs::Union{Symbol,Expr}
	rhs::Expr
end

#a general model structure with data added to the model
struct LMM
	model::Expr
	data::String
end

struct ResponseTerm
    col::Symbol
end

#make n a type with where{T} in ...
struct ConstantTerm
    n::Int64	
end

struct DataTerm <: FixedEffect
    col::Symbol 	
end

struct InteractionTerm
    cols::Vector{Symbol}	
end

struct FunctionTerm
    fun::Symbol
    cols::Union{Symbol,Expr} 	
end

struct PedigreeTerm
    var::Union{Symbol,Expr}
    path::Union{Matrix{Float64},String}
end

#runtime equivalent but does nothing
PED(var::Union{Symbol,Expr},path::Union{Matrix{Float64},String}) = PedTerm(var,path)

struct GenomicTerm
    name::Symbol
    path::String
    map::String
end

"""
        function SNP(name,path;map="")
* Defines SNP information for further analysis.
* `path` is the path for the marker file, or the matrix for the marker genotypes.
* Marker files are currently expected to be ordered as the phenotype data.
* The method to be applied to the data, `GBLUP`, `BayesB`, `BayesR` etc, is defined in the prior setting.
* Map file is optional. If not provided, a Bayesian Regression model with common variance for all SNPs will be applied. If provided, shoul match the order in the genotype file.
* One most avoid overlapping marker sets by using different `name`s.
"""
SNP(name::Symbol,path::Union{Matrix{Float64},String};map::String="") = GenomicTerm(name,path,map)
SNP(name::Symbol,path::Union{Matrix{Float64},String},map::String="") = GenomicTerm(name,path,map)
SNP(name::Symbol,path::Symbol,map::String="") = GenomicTerm(name,path,map)

struct BayesPRType
    r::Int
    v::Union{Matrix{Float64},Float64}
    name::String
    params::Bool
end

"""
        function BayesPR(r,v)
* `r` is the region size. In other words, the number of SNPs that share a common variance.
    * `1`: each SNP has its own (co)variance
    * `99`: SNPs on the same chromosome has the same (co)variance
    * `9999`: All SNPs have the same (co)variance
    * One can define any other region size, for example, 30, 40 or 100
* `v` is an estimate of the variance for the distribution of SNPs
"""
BayesPR(r::Int,v::Union{Matrix{Float64},Float64};name="BayesPR",params=false) = BayesPRType(r,v,name,params)


struct BayesBType
    pi::PiTypes
    v::VarCovarTypes
    name::String
    estimatePi::Bool
    params::Bool
end

"""
        function BayesB(pi,v;estimatePi=false)
* `pi` is the proportion of SNPs to be included in the model at each McMC cycel. If `estimatePi=true`, it is only used as a starting value.
* `v` is the variance for the prior distribution of SNPs.
* `estimatePi` is `true`if `pi` is estimated. By default it is ´false´
"""
BayesB(pi::PiTypes,v::VarCovarTypes;name="BayesB",estimatePi::Bool=false,params::Bool=false) = BayesBType(pi,v,name,estimatePi,params)

struct BayesCType
    pi::PiTypes
    v::VarCovarTypes
    name::String
    estimatePi::Bool
    params::Bool
end

"""
        function BayesC(pi,v;estimatePi=false)
* `pi` is the proportion of SNPs to be included in the model at each McMC cycel. If `estimatePi=true`, it is only used as a starting value.
* `v` is the variance for the prior distribution of SNPs.
* `estimatePi` is `true`if `pi` is estimated. By default it is ´false´
"""
BayesC(pi::PiTypes,v::VarCovarTypes;name="BayesC",estimatePi::Bool=false,params::Bool=false) = BayesCType(pi,v,name,estimatePi,params)

struct BayesRType
    pi::PiTypes
    class::Vector{Float64}
    v::VarCovarTypes
    name::String
    estimatePi::Bool
    params::Bool
end

"""
        function BayesR(pi,class,v;estimatePi=false)
* `pi` is the vector of proportion of SNPs for each variance class. If `estimatePi=true`, it is only used as a starting value.
* `class` is the vector of scales of common SPN variance for each variance class. The scales should be in the increasing order. For example, [0.0,0.0001,0.001,0.01].
* `v` is the variance for the prior distribution of SNPs.
* `estimatePi` is `true`if `pi` is estimated. By default it is ´false´
"""
BayesR(pi::PiTypes,class::Vector{Float64},v::VarCovarTypes;name="BayesR",estimatePi::Bool=false,params::Bool=false) = BayesRType(pi,class,v,name,estimatePi,params)

struct BayesRCType
    pi::PiTypes
    class::Vector{Float64}
    v::VarCovarTypes
    annot::Matrix{Int64}
    name::String
    estimatePi::Bool
    params::Bool
end

"""
        function BayesRCπ(pi,class,v,annot;estimatePi=false)
* `pi` is the vector of proportion of SNPs for each variance class. If `estimatePi=true`, it is only used as a starting value.
* `class` is the vector of scales of common SPN variance for each variance class. The scales should be in the increasing order. For example, [0.0,0.0001,0.001,0.01].
* `v` is the variance for the prior distribution of SNPs.
* `annnot` is a matrix of annottations. Rows are SNPs, columnns are annnotation classes. Entries are binary coded, 0/1.
* `estimatePi` is `true`if `pi` is estimated. By default it is ´false´
"""
BayesRCπ(pi::PiTypes,class::Vector{Float64},v::VarCovarTypes,annot::Matrix{Int64};name="BayesRCπ",estimatePi::Bool=false,params::Bool=false) = BayesRCType(pi,class,v,annot,name,estimatePi,params)
BayesRCplus(pi::PiTypes,class::Vector{Float64},v::VarCovarTypes,annot::Matrix{Int64};name="BayesRCplus",estimatePi::Bool=false,params::Bool=false) = BayesRCType(pi,class,v,annot,name,estimatePi,params)


struct BayesLogVarType
    v::Union{Matrix{Float64},Float64}
    f::LMM
    covariates::DataFrame
    varZeta::Float64
    name::String
    estimateVarZeta::Union{Float64,Bool}
end

"""
        function BayesLV(v,f,covariates,zeta)
* `v` is the variance for the prior distribution of SNPs.
* `f` is the model formula for the variance
* `covariates`is the `DataFrame` that includes explanatory varibles for the variance of each SNP.
* `zeta` is the error variance in the model for the SNP variances.
* If `estimateVarZeta` is `true`, it assumes that the error variance for the model for the SNP variances is 0.01 percent of the variance of the log-variances. If `estimateVarZeta` is `false`, it uses `zeta` as the error variance for the log-linear variances (fixed value). If `estimateVarZeta` is a `Float64`, it assumes that the error variance for the model for the SNP variances is `estimateVarZeta` percent of the variance of the log-variances.
"""
BayesLV(v::Float64,f::LMM,covariates::DataFrame,varZeta::Float64;name="BayesLV",estimateVarZeta::Union{Float64,Bool}=false) = BayesLogVarType(v,f,covariates,varZeta,name,estimateVarZeta)

struct RandomEffectType
    str::Any
    v::Union{Matrix{Float64},Vector{Float64},Float64}
    type::Int
end

"""
        function Random(str,v)
* `str` is either `I` for identity matrix, or `D` for weighted residuals
* `v` is an estimate of the variance
"""
Random(str::Any,v::Union{Matrix{Float64},Float64};type=1) = RandomEffectType(str,v,type)


struct SummaryStatistics
    m::Union{Vector{Float64},Float64}
    v::Union{Vector{Float64},Matrix{Float64},Float64}
end


