using Documenter
using NextGP

push!(LOAD_PATH,"../src/")
makedocs(
         modules  = [NextGP],
         doctest  = false,
         clean    = true,
         sitename = "NextGP.jl",
         authors  = "Emre Karaman"
         pages = Any[
            "Index" => "index.md",
            "Equations" => "equations.md",
            
         "Examples => Any[   
            "Simple example1" => "example1/example1.md",
            "Simple example2" => "example1/example1.md"
            ],
            
         ],
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information..
deploydocs(
    repo = "github.com/dataScienceToolKit/NextGP.jl.git",
    target = "build"
)
