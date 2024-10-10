using Documenter
using NextGP

push!(LOAD_PATH,"../src/")
makedocs(
         modules  = [NextGP],
         doctest  = true,
	 checkdocs=:none,
         clean    = true,
         sitename = "NextGP.jl",
#repo = GitHub("datasciencetoolkit", "NextGP.jl"),
         authors  = "Emre Karaman",
         pages = [
            "Home" => "home.md",
	    "Modules" => "index.md",
	    "Examples" => ["Example" => "Example/Example.md",
			   "PBLUP" => "PBLUP/PBLUP.md",
			   "BWGR" => "BWGR/BWGR.md",
			   "Multiple marker sets" => "MultipleMarkerSets/MultipleMarkerSets.md",
			   "BayesLV" => "BayesLV/BayesLV.md",
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
#devurl = "dev_1.0.0",
)
