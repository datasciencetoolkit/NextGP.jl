
# NextGP.jl: Modules and Functions

* `NextGP.jl` relies on several other Julia packages, which could be easily avoided and may be removed in the later releases.
* Currently one core dependency is the [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package, for model expression, and fixed effect design matrix generation.
* Checking `StatsModels.jl`'s manual for at least [`formula`](https://juliastats.org/StatsModels.jl/latest/formula/#The-@formula-language)  and  [`categorical data`](https://juliastats.org/StatsModels.jl/latest/contrasts/#Modeling-categorical-data) could be useful. 

---

## Basic Model

`NextGP.jl` uses the following basic model:


$$
\mathbf{y}=\mathbf{X}\mathbf{b}+\sum{\mathbf{Zu}}+\sum{\mathbf{M\beta}}+\mathbf{e}
$$

**b** is a vector of fixed effects 

**u** is a vector of polygenic effects in the model with a covariance
structure $\mathbf{u}\sim N\left(0,\mathbf{A}\sigma_{a}^{2}\right)$

**A** is a relationship matrix from the pedigree

**Z** is a design matrix allocating animals to records

**M** are matrices of genotypes

$\mathbf{\beta}$ are vectors  of marker effects

---

## Public Functions

Public functions of NextGP are documented here


```@docs
prep
makeA
makePed
makeG
SNP
BayesPR
```

## Internals

Functions that are used internally are documented here
