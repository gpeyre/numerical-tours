for i in 1:n
    D = min.(D, repeat(D[:,i], outer=(1,n)) + repeat(D[i,:], outer=(1,n))')
end
