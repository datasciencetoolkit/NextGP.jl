include("types.jl")

#macro model(expr::Expr)
#	m = LMM(expr,expr.args...)
#	return m
#end

#import Base.show #also export!!!
#function show(io::IO, m::LMM)
#		println("MODEL: \n $(m.model) \nLHS: \n $(m.lhs) \nRHS:")
#		for term in filter(!in(preserved), m.rhs.args)
#		    	println(" $term")
#		end
#end

macro model(expr::Expr,data::String)
      M = LMM(expr,data)
      m = lmm(M.data,M.model,M.model.args...)
      return m
end

function modelType(model::Tuple)
	models = ()
	for mi in model	
      		M = LMM(mi.model,mi.data) #expr,data
      		m = lmm(M.data,M.model,M.model.args...)
		models = (models...,m)
	end
	return models
end

import Base.show #also export!!!
function show(io::IO, m::lmm)
		println("MODEL: \n $(m.model) \nLHS: \n $(m.lhs) \nRHS:")
		for term in filter(!in(preserved), m.rhs.args)
		    	println(" $term")
		end
end

preserved = [:*,:+, :~, :-, :|, :/]
isacall(exp::Expr) = (exp.head == :call)
#isacall(exp::Expr) = ((exp == :Symbol) && !in(exp,preserved))
isacall(exp::Int) = false
isacall(exp::Symbol) = false

#multi-trait
function getLHSTerms(f;pre=preserved)
	modelLHSTerms = Dict()
	if isa(f.lhs,Symbol) 
		modelLHSTerms[f.lhs] = ResponseTerm(f.lhs)
	elseif isa(f.lhs,Expr)
		for term in f.lhs.args
			modelLHSTerms[term] = ResponseTerm(term)
		end
	else throw(DomainError("Invalid response variable"))
	end
	return modelLHSTerms
end

function getRHSTerms(f;pre=preserved)
	modelRHSTerms = Dict()
	for term in filter(!in(preserved), f.rhs.args)
		!isacall(term) && isa(term,Int) ? modelRHSTerms[Symbol("Intercept$term")] = ConstantTerm(term) : nothing
		!isacall(term) && isa(term,Symbol) ? modelRHSTerms[term] = DataTerm(term) : nothing
		isacall(term) && isdefined(Base, term.args[1]) && (getproperty(Main, term.args[1]) isa Function) && (getproperty(Main, term.args[1]) == *) ? modelRHSTerms[term] = InteractionTerm(term.args[2:end]) : nothing
		isacall(term) && isdefined(Base, term.args[1]) && (getproperty(Main, term.args[1]) isa Function) && (getproperty(Main, term.args[1]) != *) ? modelRHSTerms[term] = FunctionTerm(term.args[1],term.args[2]) : nothing
		isacall(term) && (term.args[1] == :PED) ? modelRHSTerms[term.args[2]]=PED(term.args[2:end]...) : nothing
		isacall(term) && (term.args[1] == :SNP) ? modelRHSTerms[term.args[2]]=SNP(term.args[2:end]...) : nothing	
	end
	return modelRHSTerms
end
