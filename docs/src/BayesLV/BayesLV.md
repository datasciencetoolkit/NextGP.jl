
# Bayesian log-linear variance model


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
f = @formula(y ~ 1 + SNP(M,"../data/pureGenoHOL_ref"))
```


    FormulaTerm
    Response:
      y(unknown)
    Predictors:
      1
      (M)->SNP(M, "../data/pureGenoHOL_ref")


```julia
data_LV = CSV.read(path2Data*"GWAS$(Breed)_ref",DataFrame)
```

```julia
f_LV = @formula(0 ~ x1 + x2)
```


    FormulaTerm
    Response:
      0
    Predictors:
      x1(unknown)
      x2(unknown)


```julia
priorVar = Dict(:M => BayesLV(0.001,f_LV,data_LV),
                :e => Random("I",150.0));
```


```julia
runLMEM(f,pheno,50000,10000,10;VCV=priorVar)
```



