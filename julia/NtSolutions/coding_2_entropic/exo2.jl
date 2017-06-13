e_bound = - sum(h.*log2([max(e,1e-20) for e in h]))
println("Entropy bound = $e_bound")
println("---")
Err = []

for k in 1:10
    #define constants
    m_token = 2

    #new size
    n1 = (Base.div(n,k)+1)*k

    #new vector
    x1 = zeros(n1)
    x1[1:length(x)] = x
    x1[length(x)+1:end] = 1
    x1 = x1 - 1
    x2 = []

    mult = [m_token^i for i in 0:k-1]
    for i in 1:k:n1-1
        append!(x2,sum(x1[i:i+k-1].*mult)+1)
    end

    #new probability distribution
    H = h
    for i in 1:k-1
        H = kron(H, h)
    end

    #build Huffman tree
    m = length(H)
    T=Array{Any,1}(zeros(m))

    for i in 1:m
        T[i] = (H[i],string(i))
    end

    while length(T) > 1
        sort!(T)
        t = tuple(T[1:2]...)
        q = T[1][1] + T[2][1]
        T = [T[3:end]; [(q,t)]]
    end
    T = trim(T[1])

    #find the codes
    codes = Dict()
    huffman_gencode(T,codes,"")

    #encode
    y = ""
    for e in x2
        y = string(y, codes[string(Int(e))])
    end
    #append error
    err = length(y)/length(x)
    println("Huffman(block size = $k) = $err")
    append!(Err, [err-e_bound])
end

figure(figsize = (7,5))

plot(Err, linewidth = 2)
title("Huffman block coding performance")
xlabel("Block size q")
ylabel("Huffman error - Entropy bound")
