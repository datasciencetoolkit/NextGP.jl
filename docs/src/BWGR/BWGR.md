
# Bayesian whole genome regression


```julia
using DataFrames, CSV, StatsModels, StatsBase, NextGP
```

```julia
path2Data = "../data/"
```

```julia
pheno = CSV.read(path2Data*"pheno$(Pop)_ref",DataFrame)
```


```julia
f = @formula(y ~ 1 + SNP(M,"../data/pureGenoHOL_ref","../data/map.txt"))
```


    FormulaTerm
    Response:
      y(unknown)
    Predictors:
      1
      (M)->SNP(M, "../data/pureGenoHOL_ref", "../data/map.txt")


```julia
priorVar = Dict(:M => BayesPR(9999,0.001),
                :e => Random("I",150.0));
```

```julia
runLMEM(f,pheno,50000,10000,10;VCV=priorVar)
```

    
     ---------------- Summary of input ---------------- 
    
    |Variable   Term                  Type              Levels  |
    |----------|---------------------|-----------------|--------|
    | 1        | ConstantTerm{Int64} | Vector{Float64} | 1      |
    | M        | Marker Effect              | Matrix{Float64} | 12414  |
    
    prior var-cov structure for "e" is either empty or "I" was given. An identity matrix will be used
    
     ---------------- Summary of analysis ---------------- 
    
    |Effect   Type              Str          df    scale   |
    |--------|-----------------|------------|-----|--------|
    | M      | Random (Marker) | 1 block(s) | 4.0 | 0.0005 |
    | e      | Random          | I          | 4.0 | 75.0   |


MCMC progress... 100%|███████████████████████████████████| Time: 0:10:16


