
# NextGP.jl: Modules and Functions

* `NextGP.jl` relies on several other Julia packages, which could be easily avoided and may be removed in the later releases.
* Currently one core dependency is the [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package, for model expression, and fixed effect design matrix generation.
* Checking `StatsModels.jl`'s manual for at least [`formula`](https://juliastats.org/StatsModels.jl/latest/formula/#The-@formula-language)  and  [`categorical data`](https://juliastats.org/StatsModels.jl/latest/contrasts/#Modeling-categorical-data) could be useful. 

## Basic Model

`NextGP.jl` uses the following basic model:

\begin{equation}
\mathbf{y}= \mathbf{X}\mathbf{b} + \sum_{i}\mathbf{Z}_{i}\mathbf{u}_{i}  + \sum_{j}\mathbf{M}_{j}\boldsymbol{\beta}_{j} + \mathbf{e}
\end{equation}

* $\mathbf{y}$ is a vector of phenotypes corrected
* $\mathbf{X}$ is a matrix of fixed effects
* $\mathbf{b} is a vector of fixed effects
* $\mathbf{Z_i}$ are matrices of random effects
* $\mathbf{u_i}$ are vectors of random effects
* $\mathbf{M}_{j}$ are matrices of random effects
* $\boldsymbol{\beta}_j$ are vectors of random marker effects
* $\mathbf{e}$ is the vector of random environmental
effects

## Public Functions

Public functions of NextGP are documented here


```@docs
mme
makeA
makePed
```

## Internals

Functions that are used internally are documented here
