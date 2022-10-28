
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

To install the latest official version, please use the following standard Julia command.

```@example
using Pkg
Pkg.add("NextGP")
```
## Unstable

To install the latest unofficial version (1.0.0), please use the following.

```@example
using Pkg
pkg> Pkg.add(url = "https://github.com/datasciencetoolkit/NextGP.jl", rev="dev_1.0.0")
```

