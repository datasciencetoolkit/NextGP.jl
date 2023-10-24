
# PBLUP

```julia
using CSV, StatsModels, DataFrames, NextGP
```

```julia
data = CSV.read("YOURPATH/phenotypes.csv",DataFrame)
```

```text
ID     Sire  Dam   Herds Pen BW
QGG5   QGG1  QGG2  1     1   35.0
QGG6   QGG3  QGG2  1     2   20.0
QGG7   QGG4  QGG6  1     2   25.0
QGG8   QGG3  QGG5  1     1   40.0
QGG9   QGG1  QGG6  2     1   42.0
QGG10  QGG3  QGG2  2     2   22.0
QGG11  QGG3  QGG7  2     2   35.0
QGG12  QGG8  QGG7  3     2   34.0
QGG13  QGG9  QGG2  3     1   20.0
QGG14  QGG3  QGG6  3     2   40.0
```


```julia
pedigree = "YOURPATH/pedigreeBase.txt"
```

    "YOURPATH/pedigreeBase.txt"

```text
#Pedigree for the above example
QGG1 0 0
QGG2 0 0
QGG3 0 0
QGG4 0 0
QGG5 QGG1 QGG2
QGG6 QGG3 QGG2
QGG7 QGG4 QGG6
QGG8 QGG3 QGG5
QGG9 QGG1 QGG6
QGG10 QGG3 QGG2
QGG11 QGG3 QGG7
QGG12 QGG8 QGG7
QGG13 QGG9 QGG2
QGG14 QGG3 QGG6
```

```julia
f = @formula(BW ~ Herds + Pen + PED(ID) + PED(Dam) + (1|Dam))
```


    FormulaTerm
    Response:
      BW(unknown)
    Predictors:
      Herds(unknown)
      Pen(unknown)
      (ID)->PED(ID)
      (Dam)->PED(Dam)
      (Dam)->1 | Dam



```julia
myHints = Dict(:Dam => StatsModels.FullDummyCoding(),
            :ID => StatsModels.FullDummyCoding(),
            :Herds => StatsModels.DummyCoding(),
            :Pen => StatsModels.FullDummyCoding())
```


```julia
blk = [(:Herds,:Pen)]
```


```julia
priorVar = Dict(:ID      => Random("A",150.0),
                :Dam     => Random("A",90.0),
                :(1|Dam) => Random("I",40.0),
                :e       => Random("I",350.0));
```


```julia
@time runLMEM(f,data,100000,10000,10;myHints=myHints,blockThese=blk,VCV=priorVar,userPedData=pedigree)
```

    Output folder outMCMC exists. Removing its content
    
     ---------------- Summary of input ---------------- 
    
    |Variable   Term   Type              Levels
    |----------|------|-----------------|--------|
    | Herds    | Term | Matrix{Float64} | 2      |
    | Pen      | Term | Matrix{Float64} | 2      |
    | ID       | PED  | Matrix{Bool}    | 14     |
    | Dam      | PED  | Matrix{Bool}    | 14     |
    | 1 | Dam  | |    | Matrix{Float64} | 4      |
    
    prior var-cov structure for "e" is either empty or "I" was given. An identity matrix will be used
    prior var-cov structure for ID is A. Computed A matrix (from pedigree file) will be used
    prior var-cov structure for Dam is A. Computed A matrix (from pedigree file) will be used
    prior var-cov structure for 1 | Dam is either empty or "I" was given. An identity matrix will be used
    
     ---------------- Summary of analysis ---------------- 
    
    |Effect    Type     Str   df    scale
    |---------|--------|-----|-----|-------|
    | ID      | Random | A   | 4.0 | 75.0  |
    | Dam     | Random | A   | 4.0 | 45.0  |
    | 1 | Dam | Random | I   | 4.0 | 20.0  |
    | e       | Random | I   | 4.0 | 175.0 |


    MCMC progress... 100%|███████████████████████████████████| Time: 0:00:10




