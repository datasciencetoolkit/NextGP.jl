module VCV

struct ranVar
   Str::Array{Float64, 2}
   vcov
end



function sampleRanVar(νS_ranVar,effVec,df_ranVar,Str,n)
    return((νS_ranVar + effVec'*Str*effVec)/rand(Chisq(df_ranVar + n)))
end

end



