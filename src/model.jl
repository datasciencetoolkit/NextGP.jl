include("types.jl")

macro model(expr::Expr)
	m = LMM(expr,expr.args...)
	return m
end

import Base.show #also export!!!
function show(io::IO, m::LMM)
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
function getLhsTerms(f;pre=preserved)
	modelLhsTerms = Dict()
	if is(f.lhs,Symbol) 
		modelLhsTerms[f.lhs] = ResponseTerm(f.lhs)
	elseif is(f.lhs,Expr)
		for term in f.lhs.args
			modelLhsTerms[term] = ResponseTerm(term)
		end
	else throw(DomainError("Invalid response variable"))
	end
end

function getRhsTerms(f;pre=preserved)
	modelRhsTerms = Dict()
	for term in filter(!in(preserved), f.rhs.args)
		!isacall(term) && (term==1) ? modelRhsTerms[:(Intercept)] = ConstantTerm(term) : nothing
		!isacall(term) && isa(term,Symbol) ? modelRhsTerms[term] = DataTerm(term) : nothing
		isacall(term) && isdefined(Base, term.args[1]) && (getproperty(Main, term.args[1]) isa Function) && (getproperty(Main, term.args[1]) == *) ? modelRhsTerms[term] = InteractionTerm(term.args[2:end]) : nothing
		isacall(term) && isdefined(Base, term.args[1]) && (getproperty(Main, term.args[1]) isa Function) && (getproperty(Main, term.args[1]) != *) ? modelRhsTerms[term] = FunctionTerm(term.args[1],term.args[2]) : nothing
		isacall(term) && (term.args[1] == :PED) ? modelRhsTerms[term.args[2]]=PED(term.args[2:end]...) : nothing
		isacall(term) && (term.args[1] == :SNP) ? modelRhsTerms[term.args[2]]=SNP(term.args[2:end]...) : nothing	
	end
	return modelRhsTerms
end
