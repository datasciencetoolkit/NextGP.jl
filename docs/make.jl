using Documenter
using NextGP

makedocs(
         root     = "./",
	 source   = "src",
         modules  = [NextGP],
         doctest  = false,
         clean    = true,
         sitename = "NextGP.jl",
         authors  = "Emre Karaman",
         pages = [
            "Index" => "index.md",
         ],
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information..
deploydocs(
    repo = "github.com/datasciencetoolkit/NextGP.jl.git",
    branch = "main",
)
