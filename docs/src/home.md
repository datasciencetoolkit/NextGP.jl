
# About

## What is available

* Currently, only univariate analysis are implemented
* Rely on [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package for model formulation and fixed effects definitions
* Any kind of random effects including random marker effects                                 
  - Additive genetic effects
  - Maternal effects
  - Permanent environmental effects
* Bayesian whole-genome regression methods
  - BayesPR (BayesA<->BRR)
  - BayesB
  - BayesC
  - BayesR
  - BayesLV
* Correlated marker effects (only for BayesPR)
* Unknown (co)variance components (e.g., marker,additive genetic,residual...)
* It is users responsibility to make sure that the order of the individuals in the data sets aligns

## Pkg Registry

To install the latest official version, please use the following standard Julia command.

```@example
using Pkg
Pkg.add("NextGP")
```
## Unstable

To install the latest unofficial version, please use the following.

```@example
using Pkg
Pkg.add(url = "https://github.com/datasciencetoolkit/NextGP.jl", rev="dev_1.2.0")
```

