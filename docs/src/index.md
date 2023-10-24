
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
runLMEM
makeA
makePed
makeG
SNP
BayesPR
BayesB
BayesC
BayesR
BayesLV
```

## Internals

Functions that are used internally are documented here

```@docs
prep
```

## Convenience functions

### Convergency checking

```julia
using MCMCChains,StatsPlots,StatsBase,CSV

function summaryMCMC(param;summary=false,plots=false,outFolder=pwd()*"/outMCMC")
        param = CSV.read("$outFolder/$(param)Out",DataFrame,header=true)
        namesParam = names(param)
        param = Matrix(param)
                if summary==true
                        chn = Chains(param,namesParam)
                        display(chn)
                        if plots==true
                                display(plot(chn))
                        end
                        param = mean(Matrix(param),dims=1)
                else
                        param = mean(Matrix(param),dims=1)
                end
        return param
end

```

* If `summary=true`, will print convergency statistics for McMC
* If `plots=true`, will print trace plot(s) of McMC
* `outFolder` is the folder for the McMC output. By default it looks for the folder "outMCMC" in the current directory.

---


