
# Multiple marker sets


```julia
Breed = "HOL"
```

```julia
using DataFrames, CSV, StatsModels, StatsBase, NextGP
```


```julia
myHints = Dict(:lact => StatsModels.FullDummyCoding(),
               :herd => StatsModels.FullDummyCoding())
```

    Dict{Symbol, StatsModels.FullDummyCoding} with 3 entries:
      :herd => FullDummyCoding()
      :lact => FullDummyCoding()

```julia
f = @formula(y ~ 1 + lact + herd + dim + wilmink + SNP(A,"arc.txt") + SNP(B,"bac.txt"))
```

    FormulaTerm
    Response:
      y(unknown)
    Predictors:
      1
      lact(unknown)
      herd(unknown)
      dim(unknown)
      wilmink(unknown)
      (A)->SNP(A, "arc.txt")
      (B)->SNP(B, "bac.txt")




```julia
blk = [(Symbol(1),:lact,:herd)]
```

    1-element Vector{Tuple{Symbol, Symbol}}:
     (Symbol("1"), :lact)


```julia
priorVar = Dict(:A => BayesPR(9999,0.04),
                :B => BayesPR(9999,0.04),
                :e => Random("I",2500.0));
```


```julia
runLMEM(f,pheno,5000,500,10;outFolder="pure$Breed",VCV=priorVar,myHints=myHints,blockThese=blk)
```

    Output folder pureHOL exists. Removing its content
    
     ---------------- Summary of input ---------------- 
    
    |[1m Variable [0m|[1m Term                [0m|[1m Type            [0m|[1m Levels [0m|
    |----------|---------------------|-----------------|--------|
    | 1        | ConstantTerm{Int64} | Vector{Float64} | 1      |
    | lact     | Term                | Matrix{Float64} | 6      |
    | herd     | Term                | Matrix{Float64} | 6      |
    | dim      | Term                | Vector{Float64} | 1      |
    | wilmink  | Term                | Vector{Float64} | 1      |
    | A        | Marker Effect       | Matrix{Float64} | 189    |
    | B        | Marker Effect       | Matrix{Float64} | 3894   |
    [32mprior var-cov structure for "e" is either empty or "I" was given. An identity matrix will be used[39m
    [32mNo map was provided. Running Bayesian Random Regression (BRR) with all SNP as 1 region[39m
    [32mNo map was provided. Running Bayesian Random Regression (BRR) with all SNP as 1 region[39m
    
     ---------------- Summary of analysis ---------------- 
    
    |[1m Effect [0m|[1m Type            [0m|[1m Str        [0m|[1m df  [0m|[1m scale  [0m|
    |--------|-----------------|------------|-----|--------|
    | A      | Random (Marker) | 1 block(s) | 4.0 | 0.02   |
    | B      | Random (Marker) | 1 block(s) | 4.0 | 0.02   |
    | e      | Random          | I          | 4.0 | 1250.0 |


