e_bound = - np.sum(h*np.log2([max(e,1e-20) for e in h]))
print("Entropy bound = %f" %e_bound)
print("---")
Err = []

for k in range(1,11):
    #define constants
    m_token = 2
    
    #new size
    n1 = (n//k+1)*k
    
    #new vector
    x1 = np.zeros(n1)
    x1[:len(x)] = x
    x1[len(x):] = 1
    x1 = x1 - 1
    x2 = []
    
    for i in range(0,n1,k):
        mult = [m_token**i for i in range(k)]
        x2.append(sum(x1[i:i+k]*mult))
        
    #new probability distribution
    H = h
    for i in range(1,k):
        H = np.kron(H, h)

    #build Huffman tree
    m = len(H)
    T = [0] * m 
    for i in range(m):
        T[i] = (H[i],str(i))
    
    while len(T) > 1: 
        T.sort()
        t = tuple(T[:2])
        q = T[0][0] + T[1][0]
        T = T[2:] + [(q,t)]
        
    T = trim(T[0])

    #find the codes
    codes = {}
    huffman_gencode(T,codes,"") 
    
    #encode
    y = ""
    for e in x2:
        y = y + codes[str(int(e))]
        
    #append error
    err = len(y)/len(x)
    print("Huffman(block size = %i) = %f" %(k,err))
    Err.append(err-e_bound)

plt.figure(figsize = (7,5))

plt.plot(Err, linewidth = 2)
plt.title("Huffman block coding performance")
plt.xlabel("Block size q")
plt.ylabel("Huffman error $-$ Entropy bound")

plt.show()