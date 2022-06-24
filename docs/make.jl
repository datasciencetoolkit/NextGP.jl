using Documenter
using NextGP

push!(LOAD_PATH,"../src/")

makedocs(
    sitename = "NextGP",
    format = Documenter.HTML(prettyurls = false),
    modules = [NextGP]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
