using Documenter
using NextGP

push!(LOAD_PATH,"../src/")
makedocs(
         modules  = [NextGP],
         doctest  = true,
         clean    = true,
         sitename = "NextGP.jl",
         authors  = "Emre Karaman",
         pages = [
            "Home" => "home.md",
	    "Modules" => "index.md",
	    "Examples" => ["Example1" => "example1/example1.md",
			   "Example2" => "example2/example2.md",
			],
	    "Citation" => "citation.md",
         ],
	format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true"
    )
)
# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information..
deploydocs(
    repo = "github.com/datasciencetoolkit/NextGP.jl",
)
