using Documenter
using NextGP

push!(LOAD_PATH,"../src/")
makedocs(sitename="NextGP.jl Documentation",
         pages = [
            "Index" => "index.md",
            "Equations" => "equations.md",
         ],
         format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true")
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/dataScienceToolKit/NextGP.jl.git",
    branch = "gh-pages"
)
