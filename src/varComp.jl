#set up (co)variance structures for E
function setVarCovStrE!(eSet::Symbol,E,priorVCV,nData,varE)	
	#no inverse implemented yet!
	if haskey(priorVCV,eSet)	
		if isempty(priorVCV[eSet].str) || priorVCV[eSet].str=="I" 
				printstyled("prior var-cov structure for \"e\" is either empty or \"I\" was given. An identity matrix will be used\n"; color = :green)
				E[eSet][:str] = "I"
				E[eSet][:iVarStr] = [] #Matrix(1.0I,nData,nData)
				priorVCV[eSet] = Random("I",priorVCV[eSet].v)
		elseif isa(priorVCV[eSet].str,Vector) # D
				E[eSet][:str] = "D"
				E[eSet][:iVarStr] = inv.(priorVCV[eSet].str) #inv(Diagonal(priorVCV[eSet].str))
#				error("var-cov structure \"D\" has not been implemented yet")
				printstyled("prior var-cov structure for \"e\" is \"D\". User provided \"D\" matrix (d_ii = 1/w_ii) will be used\n"; color = :green)
		else 
				error("provide a valid prior var-cov structure (\"I\", \"D\" or leave it empty \"[]\") for \"e\" ")
		end
	else	
		printstyled("prior var-cov for \"e\" is fully  empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		E[eSet][:str] = "I"
		E[eSet][:iVarStr] = [] #Matrix(1.0I,nData,nData)
		#just add to priors
		priorVCV[eSet] = Random("I",100.0)
	end
	varE[eSet] = priorVCV[eSet].v
end

function setVarCovStrE!(eSet::Tuple{Vararg{Symbol}},E,priorVCV,nData,varE)	
	#no inverse implemented yet!
	if haskey(priorVCV,eSet)	
		if isempty(priorVCV[eSet].str) || priorVCV[eSet].str=="I" 
				printstyled("prior var-cov structure for \"e\" is either empty or \"I\" was given. R0 (kron) I will be used\n"; color = :green)
				E[eSet][:str] = "I"
				E[eSet][:iVarStr] = [] #Matrix(1.0I,nData,nData)
				priorVCV[eSet] = Random("I",priorVCV[eSet].v)
		elseif isa(priorVCV[eSet].str,Vector) # D
				E[eSet][:str] = "D"
				E[eSet][:iVarStr] = inv.(priorVCV[eSet].str) #inv(Diagonal(priorVCV[eSet].str))
#				error("var-cov structure \"D\" has not been implemented yet")
				printstyled("prior var-cov structure for \"e\" is \"D\". User provided \"D\" matrix (d_ii = 1/w_ii) will be used\n"; color = :green)
		else 
				error("provide a valid prior var-cov structure (\"I\", \"D\" or leave it empty \"[]\") for \"e\" ")
		end
	else	
		printstyled("prior var-cov for \"e\" is fully  empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		E[eSet][:str] = "I"
		E[eSet][:iVarStr] = [] #Matrix(1.0I,nData,nData)
		#just add to priors
		priorVCV[eSet] = Random("I",100.0)
	end
	varE[eSet] = priorVCV[eSet].v
end


###############################
#df, shape, scale...															
function varCovE!(E,priorVCV)
	for eSet ∈ keys(E)
		#df
		E[eSet][:df] = 3.0+size(priorVCV[eSet].v,1)
		#scale
                nEComp = size(priorVCV[eSet].v,1)
		#priorVCV[eSet].v is a temporary solution
		nEComp > 1 ? E[eSet][:scale] = priorVCV[eSet].v .* (E[eSet][:df]-nEComp-1.0)  : E[eSet][:scale] = priorVCV[eSet].v * (E[eSet][:df]-2.0)/E[eSet][:df] #I make float and array of float														
        end
#	if priorVCV[:e].v==0.0
#		priorVCV[:e].v  = 0.0005
#       		E[:scale]     = 0.0005
#        else
#       		E[:scale]    = priorVCV[:e].v*(E[:df]-2.0)/E[:df]    
#   	end
end

#set up (co)variance structures for U

function setVarCovStrU!(zSet::ExprOrSymbolOrTuple,Z::Dict,priorVCV,varU::Dict)
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
		varU[zSet] = priorVCV[zSet].v
	else	
		printstyled("prior var-cov for $zSet is empty. An identity matrix will be used with mean=0 and variance=100\n"; color = :green)
		varU[zSet] = 100
		priorVCV[zSet] = Random("I",100)
		Z[zSet][:iVarStr] = Matrix(1.0I,Z[zSet][:dims][2],Z[zSet][:dims][2])
	end
end

#df, shape, scale...															
function varCovZ!(Z,priorVCV)
	for zSet ∈ keys(Z)
		#df
		Z[zSet][:df] = 3.0+size(priorVCV[zSet].v,1)
		#scale
                nZComp = size(priorVCV[zSet].v,1)
		#priorVCV[zSet].v is a temporary solution
		nZComp > 1 ? Z[zSet][:scale] = priorVCV[zSet].v .* (Z[zSet][:df]-nZComp-1.0)  : Z[zSet][:scale] = priorVCV[zSet].v * (Z[zSet][:df]-2.0)/Z[zSet][:df] #I make float and array of float														
        end
end


#set up (co)variance priors for markers
#structure comes from method chosen, so it is not here.

function varCovM!(M,priorVCV,varBeta)
	#df
	for mSet ∈ keys(M)
		M[mSet][:df] = haskey(priorVCV,mSet) ? 3.0+size(priorVCV[mSet].v,1) : 4.0
	end
	#scale
	#when hyper-parameters are estimated, all models will have only one estimate. Not per SNP, for now!
	#But in general, all models get same prior anyways.
        for mSet ∈ keys(M)
		if haskey(priorVCV,mSet)
                	nMComp = size(priorVCV[mSet].v,1)
			if priorVCV[mSet].params==false
				println("hyperparameters are fixed")
				M[mSet][:params] = false
                		M[mSet][:scale] = nMComp>1 ? priorVCV[mSet].v .* (M[mSet][:df]-nMComp-1.0)  : priorVCV[mSet].v * (M[mSet][:df]-2.0)/(M[mSet][:df]) #I make array of float and float 
			elseif priorVCV[mSet].params==true
				println("hyperparameters are estimated")
				M[mSet][:params] = true
				M[mSet][:scale] = nMComp>1 ? [priorVCV[mSet].v .* (M[mSet][:df]-nMComp-1.0)]  : [priorVCV[mSet].v * (M[mSet][:df]-2.0)/(M[mSet][:df])] #I make array of float and float
			end
		else
			M[mSet][:params] = false
			nMComp = 1
			M[mSet][:scale] = 0.05 * (M[mSet][:df]-2.0)/(M[mSet][:df]) #I make float and array of float
		end

		#for storage
        	for mSet ∈ keys(M)
			if haskey(priorVCV,mSet)
                		varBeta[mSet] = [priorVCV[mSet].v for i in 1:M[mSet][:nVarCov]]
			else
			varBeta[mSet] = [0.05 for i in 1:M[mSet][:nVarCov]]
			end
        	end
        end
end




