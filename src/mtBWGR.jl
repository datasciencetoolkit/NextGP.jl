vec2Mat(d) = reshape(d,1,length(d))

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
			#M[mSet][:scale]   = []
		elseif priorVCV[mSet].name == "BayesB"
			#M[mSet][:logPi]       = log.(priorVCV[mSet].pi)
			M[mSet][:gammaComb] = collect.(Iterators.product(fill(0:1, length(mSet))...) |> collect |> vec)
			
			for g in 1:length(M[mSet][:gammaComb])
    				println("gamma $(M[mSet][:gammaComb][g]) has prior $(priorVCV[mSet].pi[g]) piDelta \n")
			end

			deltaComb = Vector(Diagonal.(collect.(M[mSet][:gammaComb]))) #this is in matrix form [1 0 0;0 1 0;0 0 1]
			println("deltaComb: $deltaComb")
			M[mSet][:deltaComb] = deltaComb

			M[mSet][:gammaComb] = vec2Mat.(M[mSet][:gammaComb])
			for g in 1:length(M[mSet][:gammaComb])
    				println("gamma $(M[mSet][:gammaComb][g]) has prior $(priorVCV[mSet].pi[g]) piDelta \n")
			end
			
			M[mSet][:method]      = "BayesB"
			M[mSet][:funct]       = sampleBayesB!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = length(theseRegions)
			#gamma and delta estimates      
			M[mSet][:gammaHat] = fill(ones(1,length(mSet)),M[mSet][:dims][2])
			M[mSet][:deltaHat]    = fill(Matrix(Diagonal(ones(length(mSet)))),M[mSet][:dims][2]) #priorVCV[mSet].pi ##this is fixed if estPi=false, or just a starting value if estPi=true
			#piDelta estimate store
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = priorVCV[mSet].pi #add when not available
			M[mSet][:vClass]      = [0 1] #2 variance class, one with own, one with null
			#M[mSet][:scale]      = []
		elseif priorVCV[mSet].name == "BayesC"
			#M[mSet][:logPi]       = log.(priorVCV[mSet].pi)
			M[mSet][:gammaComb] = collect.(Iterators.product(fill(0:1, length(mSet))...) |> collect |> vec)
			
			for g in 1:length(M[mSet][:gammaComb])
    				println("gamma $(M[mSet][:gammaComb][g]) has prior $(priorVCV[mSet].pi[g]) piDelta \n")
			end

			deltaComb = Vector(Diagonal.(collect.(M[mSet][:gammaComb]))) #this is in matrix form [1 0 0;0 1 0;0 0 1]
			println("deltaComb: $deltaComb")
			M[mSet][:deltaComb] = deltaComb

			M[mSet][:gammaComb] = vec2Mat.(M[mSet][:gammaComb])
			for g in 1:length(M[mSet][:gammaComb])
    				println("gamma $(M[mSet][:gammaComb][g]) has prior $(priorVCV[mSet].pi[g]) piDelta \n")
			end
			
			M[mSet][:method]      = "BayesC"
			M[mSet][:funct]       = sampleBayesC!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = 1
			#gamma and delta estimates      
			M[mSet][:gammaHat]    = fill(ones(1,length(mSet)),M[mSet][:dims][2])
			M[mSet][:deltaHat]    = fill(Matrix(Diagonal(ones(length(mSet)))),M[mSet][:dims][2]) #priorVCV[mSet].pi ##this is fixed if estPi=false, or just a starting value if estPi=true
			#piDelta estimate store
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = priorVCV[mSet].pi #add when not available
			M[mSet][:vClass]      = [0 1] #2 variance class, one with own, one with null
			#M[mSet][:scale]      = []
		end
	end
end
