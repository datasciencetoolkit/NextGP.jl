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

function getTerms(f;pre=preserved)
	modelTerms = Dict()
	for term in filter(!in(preserved), f.rhs.args)
		!isacall(term) && (term==1) ? modelTerms[:(Intercept)] = ConstantTerm(term) : nothing
		!isacall(term) && isa(term,Symbol) ? modelTerms[term] = DataTerm(term) : nothing
		isacall(term) && isdefined(Base, term.args[1]) && (getproperty(Main, term.args[1]) isa Function) && (getproperty(Main, term.args[1]) == *) ? modelTerms[term] = InteractionTerm(term.args[2:end]) : nothing
		isacall(term) && isdefined(Base, term.args[1]) && (getproperty(Main, term.args[1]) isa Function) && (getproperty(Main, term.args[1]) != *) ? modelTerms[term] = FunctionTerm(term.args[1],term.args[2]) : nothing
		isacall(term) && (term.args[1] == :PED) ? modelTerms[term.args[2]]=PedTerm(term.args[2:end]...) : nothing
		isacall(term) && (term.args[1] == :SNP) ? modelTerms[term.args[2]]=GenomicTerm(term.args[2:end]...) : nothing	
	end
	return modelTerms
end
