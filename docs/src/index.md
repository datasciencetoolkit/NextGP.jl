
# NextGP.jl: Modules and Functions

* `NextGP.jl` relies on several other Julia packages, which could be easily avoided and may be removed in the later releases.
* Currently one core dependency is the [`StatsModels.jl`](https://juliastats.org/StatsModels.jl/latest/) package, for model expression, and fixed effect design matrix generation.
* Checking `StatsModels.jl`'s manual for at least [`formula`](https://juliastats.org/StatsModels.jl/latest/formula/#The-@formula-language)  and  [`categorical data`](https://juliastats.org/StatsModels.jl/latest/contrasts/#Modeling-categorical-data) could be useful. 


## Public Functions

Public functions of NextGP are documented here


```@docs
mme
makeA
makePed
```

## Internals

Functions that are used internally are documented here
