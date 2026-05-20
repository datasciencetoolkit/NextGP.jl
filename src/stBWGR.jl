function stBWGR!(M,mSet,priorVCV,beta)
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
			M[mSet][:logPi]       = [log(1.0 .- priorVCV[mSet].pi) log(priorVCV[mSet].pi)] #not fitted, fitted
#					M[mSet][:logPiIn]     = log(priorVCV[mSet].pi)
#					M[mSet][:logPiOut]    = log(1.0 .- priorVCV[mSet].pi)
			M[mSet][:method]      = "BayesB"
			M[mSet][:funct]       = sampleBayesB!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = length(theseRegions)
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = [1.0 .- priorVCV[mSet].pi priorVCV[mSet].pi] #not fitted, fitted
			M[mSet][:vClass]      = [0 1] #2 variance class, one with own, one with null
		elseif priorVCV[mSet].name == "BayesC"
			M[mSet][:logPi]       = [log(1.0 .- priorVCV[mSet].pi) log(priorVCV[mSet].pi)] #not fitted, fitted
#					M[mSet][:logPiIn]     = log(priorVCV[mSet].pi)
#					M[mSet][:logPiOut]    = log(1.0 .- priorVCV[mSet].pi)
			M[mSet][:method]      = "BayesC"
			M[mSet][:funct] = sampleBayesC!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = 1
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = [1.0 .- priorVCV[mSet].pi priorVCV[mSet].pi] #not fitted, fitted
			M[mSet][:vClass]      = [0 1] #2 variance class, one with common, one with null
		elseif priorVCV[mSet].name == "BayesR"
			M[mSet][:logPi]       = log.(priorVCV[mSet].pi)
			M[mSet][:vClass]      = priorVCV[mSet].class
			M[mSet][:method]      = "BayesR"
			M[mSet][:funct]       = sampleBayesR!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = 1
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = deepcopy(priorVCV[mSet].pi)
		elseif priorVCV[mSet].name == "BayesRCπ"
			M[mSet][:vClass]      = priorVCV[mSet].class
			M[mSet][:method]      = "BayesRCπ"
			M[mSet][:funct]       = sampleBayesRCπ!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = size(priorVCV[mSet].annot,2)
			M[mSet][:logPi]       = [log.(priorVCV[mSet].pi) for i in 1:M[mSet][:nVarCov]]
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = [priorVCV[mSet].pi for i in 1:M[mSet][:nVarCov]]
			M[mSet][:annotInput]  = deepcopy(priorVCV[mSet].annot)
			M[mSet][:annotProb]   = priorVCV[mSet].annot./sum(priorVCV[mSet].annot,dims=2)
					#If all annotations are zero, prob is NA. I make it "0" here
					#But if all zero, those SNPs should be in a seperate marker set.
					#So i cancel this
					#M[mSet][:annotProb][findall([all(iszero, row) for row in eachrow(priorVCV[mSet].annot)]),:] .= 0.0
					#
			M[mSet][:annotNonZeroPos]   = [findall(!iszero, row) for row in eachrow(priorVCV[mSet].annot)]
#					M[mSet][:annotNonZero]= getindex.(Ref(priorVCV[mSet].annot),M[mSet][:annotNonZeroPos])
			M[mSet][:annotCat]    = zeros(Int64,1,M[mSet][:dims][2])
		elseif priorVCV[mSet].name == "BayesRCplus"
			M[mSet][:vClass]      = priorVCV[mSet].class
			M[mSet][:method]      = "BayesRCplus"
			M[mSet][:funct]       = sampleBayesRCplus!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = size(priorVCV[mSet].annot,2)
			M[mSet][:logPi]       = [log.(priorVCV[mSet].pi) for i in 1:M[mSet][:nVarCov]]
			M[mSet][:estPi]       = priorVCV[mSet].estimatePi
			M[mSet][:piHat]       = [priorVCV[mSet].pi for i in 1:M[mSet][:nVarCov]]
			M[mSet][:annotInput]  = deepcopy(priorVCV[mSet].annot)
			M[mSet][:annotProb]   = priorVCV[mSet].annot./sum(priorVCV[mSet].annot,dims=2)
			M[mSet][:annotNonZeroPos]   = [findall(!iszero, row) for row in eachrow(priorVCV[mSet].annot)]
			M[mSet][:annotCat]    = zeros(Int64,1,M[mSet][:dims][2])
		elseif priorVCV[mSet].name == "BayesLV"
			M[mSet][:method]      = "BayesLV"
			M[mSet][:funct]       = sampleBayesLV!
			theseRegions          = [r:r for r in 1:M[mSet][:dims][2]]
			M[mSet][:regionArray] = theseRegions
			M[mSet][:nVarCov]     = length(theseRegions)
					#logVar can be created in a smarter way, maybe together with var??...
			M[mSet][:logVar]      = [log(priorVCV[mSet].v) for i in 1:M[mSet][:nVarCov]]

			varRHSTerms = getRHSTerms(priorVCV[mSet].f)
			varLHSTerms = getLHSTerms(priorVCV[mSet].f)
			
			println("varRHSTerms: $varRHSTerms")
			println("varLHSTerms: $varLHSTerms")
			
			if isa(priorVCV[mSet].f.data,String)
				println("reading data from: $(priorVCV[mSet].f.data)")
				varData = CSV.read(priorVCV[mSet].f.data,DataFrames.DataFrame,header=true,delim=' ',pool=false,stringtype=String)
			elseif isa(priorVCV[mSet].f.data,Symbol)
				println("using $(priorVCV[mSet].f.data)")
				varData = getproperty(Main,priorVCV[mSet].f.data)
			else throw(ArgumentError("Could not understand your data. Provide a file path or DataFrame"))
			end
			println("size of marker variance data set read: $(size(varData))")
			println("variables in the marker variance data set: $(names(varData))")

			#dMatrix = hcat([varData[!,k] for (k,v) in varRHSTerms]...)
			dMatrix = hcat([designMat(k,v,varData)[:data] for (k,v) in varRHSTerms]...)
			println("dMatrix: $dMatrix")
			
			#dMatrix               = designMat(k,v,varData)
			M[mSet][:covariates] = dMatrix
			M[mSet][:covariatesT]= transpose(dMatrix)
			M[mSet][:c]          = rand(size(dMatrix,2))
			M[mSet][:SNPVARRESID] = rand(size(dMatrix,1))
			M[mSet][:iCpC] = M[mSet][:covariatesT]*M[mSet][:covariates] #taking inverse after 
			if isa(M[mSet][:iCpC],Matrix{Float64}) 
				M[mSet][:iCpC] += Matrix(I*minimum(abs.(diag(M[mSet][:iCpC])./10000)),size(M[mSet][:iCpC]))
			end
 		    M[mSet][:iCpC]  = inv(M[mSet][:iCpC])
			M[mSet][:varZeta]  = [priorVCV[mSet].varZeta]
			M[mSet][:estVarZeta] = priorVCV[mSet].estimateVarZeta
			dMatrix = 0
		end
	end
end
