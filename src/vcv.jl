module VCV

struct ranVar
   Str::Array{Float64, 2}
   vcov
end

end


function sampleRanVar(νS_ranVar,effVec,df_ranVar,Str,n)
    return((νS_ranVar + effVec'*Str*effVec)/rand(Chisq(df_ranVar + n)))
end


function sampleEpsi!(Ai11,nY1,n1,Z11p,zpz,varE,varG,ycorr,ϵ,epsRows)
    λ = varE/varG
@views    ycorr[1:nY1] .+= ϵ[epsRows] #Z11*ϵ
    Yi = Z11p*ycorr[1:nY1]
    for i in 1:n1
        ϵ[i] = 0.0
        rhsEpsi = Yi[i] - λ*dot(view(Ai11,:,i),ϵ)
        lhsEpsi = zpz[i] + view(Ai11,i,i)*λ
        invLhsEpsi = 1.0/lhsEpsi
        meanEpsi = invLhsEpsi*rhsEpsi
        ϵ[i] = rand(Normal(meanEpsi,sqrt(invLhsEpsi*varE)))
    end
@views    ycorr[1:nY1] .-= ϵ[epsRows] #Z11*ϵ
end


