#set up (co)variance structures
function setVarCovStr!(zSet::ExprOrSymbolOrTuple,Z::Dict,priorVCV,varU_prior::Dict)
	if haskey(priorVCV,zSet)	
		if ismissing(priorVCV[zSet].str) || priorVCV[zSet].str=="I" 
			printstyled("prior var-cov structure for $zSet is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
			Z[zSet][:iVarStr] = Matrix(1.0I,Z[zSet][:dims][2],Z[zSet][:dims][2])
		elseif priorVCV[zSet].str=="A"
			printstyled("prior var-cov structure for $zSet is A. Computed A matrix (from pedigree file) will be used\n"; color = :green)
			isa(zSet,Tuple) ? Z[zSet][:iVarStr] = Z[zSet[1]][:iVarStr] : Z[zSet][:iVarStr] = Z[zSet][:iVarStr]
		elseif priorVCV[zSet].str=="G"
                        printstyled("prior var-cov structure for $zSet is G. Computed G matrix will be used\n"; color = :green)
			isa(zSet,Tuple) ? Z[zSet][:iVarStr] = Z[zSet[1]][:iVarStr] : Z[zSet][:iVarStr] = Z[zSet][:iVarStr]
		else 	Z[zSet][:iVarStr] = inv(priorVCV[zSet].str)
		end
		varU_prior[zSet] = priorVCV[zSet].v
	else	
		printstyled("prior var-cov for $zSet is empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		varU_prior[zSet] = 100
		priorVCV[zSet] = Random("I",100)
		Z[zSet][:iVarStr] = Matrix(1.0I,Z[zSet][:dims][2],Z[zSet][:dims][2])
	end
end



