# adapted from http://morotalab.org/Mrode2005/relmat/createA.txt
function makeA(s::Any, d::Any)
    s = convert(Vector{Int64},s)
    d = convert(Vector{Int64},d)
    n = length(s)
    N = n + 1
    A = zeros(N, N)
    s = (s .== 0)*N + s
    d = (d .== 0)*N + d
for i in 1:n
    A[i,i] = 1.0 + A[s[i], d[i]]/2.0
        for j in (i+1):n
            if j > n break end
                A[i,j] = ( A[i, s[j]] + A[i,d[j]] )/2.0
                A[j,i] = A[i,j]
    end
    end
return(A[1:n, 1:n])
end

