using Documenter
using NextGP

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "NextGP",
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true"),
    modules = [NextGP]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.

deploydocs(
    repo = "https://github.com/dataScienceToolKit/NextGP.jl"
)
