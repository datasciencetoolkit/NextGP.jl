
# About

## What is available

* Currently, only univariate analysis are implemented
* Rely on StatsModels.jl package for model formulation and fixed effects definitions
* Any kind of random effects including random marker effects                                 
  - Additive genetic effects
  - Maternal effects
  - Permanent environmental effects
* Correlated marker effects for multi-breed analysis
* Unknown (co)variance components (e.g., marker,additive genetic,residual...)

## Pkg Registry
```@repl
using Pkg
Pkg.add("NextGP")
```
## Unstable
```@example
pkg> Pkg.add(url = "https://github.com/datasciencetoolkit/NextGP.jl", rev="dev")
```

```@docs
MyFunction
```
