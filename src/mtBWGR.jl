function mtBWGR!(M,mSet,priorVCV,beta)
	if !haskey(priorVCV,mSet)
		M[mSet][:method] = "BayesPR"
		M[mSet][:funct] = sampleBayesPR!
		theseRegions = [1:r for r in M[mSet][:dims][2]]
		M[mSet][:regionArray] = theseRegions
		M[mSet][:nVarCov] = length(theseRegions)
	else
		if priorVCV[mSet].name == "BayesPR"
			M[mSet][:method] = "BayesPR"
			M[mSet][:funct] = sampleBayesPR!
			if isempty(M[mSet][:map])		
				if priorVCV[mSet].r == 1
					printstyled("No map was provided. Running Bayesian Random Regression (BRR) with 1 SNP region size\n"; color = :green)
					theseRegions = [r:r for r in 1:M[mSet][:dims][2]]
					M[mSet][:regionArray] = theseRegions
				elseif priorVCV[mSet].r == 9999
					printstyled("No map was provided. Running Bayesian Random Regression (BRR) with all SNP as 1 region\n"; color = :green)
					theseRegions = [1:r for r in M[mSet][:dims][2]]
					M[mSet][:regionArray] = theseRegions
				else throw(ArgumentError("Please enter a valid region size (1 or 9999) or provide a map file"))
				end
			else
				theseRegions = prep2RegionData(outPut,mSet,M[mSet][:map],priorVCV[mSet].r)
				M[mSet][:regionArray] = theseRegions
			end
			M[mSet][:nVarCov] = length(theseRegions)
			M[mSet][:scale]   = []
		end
	end
end
