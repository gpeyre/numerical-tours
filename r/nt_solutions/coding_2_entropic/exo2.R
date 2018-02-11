options(repr.plot.width=4, repr.plot.height=4)


e_bound = -sum(h*log2(c(max(h[1], 1e-20), max(h[2], 1e-20))))
print(paste("Entropy Bound =", e_bound))
print("---")
Err = c()

for (k in (2:11))
{
  #define constants
  m_token = 2
  #new size
  n1 = (floor(n/k)+1)*k
  
  #new vector
  x1 = matrix(0, nrow=1, ncol=n1)
  x1[1:length(x)] = x
  x1[length(x):length(x1)] = 1
  x1 = x1 - 1
  x2 = c()
  for (i in seq(1,n1,by=k))
  {
    mult = m**c(1:k-1)
    x2 = c(x2, sum(x1[i:(i+k)]*mult))
  }
  
  #new probability distribution
  H = h
  for (i in (1:(k-1)))
  {
    H = kronecker(H,h)
  }
  
  #build Huffman tree
  m = length(H)
  T = list(list())
  for (i in (1:m))
  {
    T[[i]] = list(toString(i), H[i])
  }
  while (length(T)>=2)
  {
    T = T[order(sapply(T,'[[',2))]
    q =  as.numeric(T[[1]][2])+as.numeric(T[[2]][2])
    t = T[1:2]
    T = T[-(1:2)]
    T[[length(T)+1]] = list(t,q)
  }
  K = list()
  K[[1]] = list(trim(T)[[1]][[1]][1],trim(T)[[2]])
  T = K 
  
  
  #find the codes
  codes = list()
  c = compute_huffcode(as.numeric(H))
  for (i in (1:length(h)))
  {
    codes[[toString(i)]] = c[i]
  }
  
  #encode
  y = ""
  
  for (e in x2)
  {
    y = paste0(y,codes[[toString(e)]])
  }
  #append error
  err = nchar(y)/length(x)
  print(paste("Huffman(block size = ", k, "=", err))
  Err = c(Err,err-e_bound)
  
} 

plot(Err, xlab="Block size q", ylab ="Huffman error $-$ Entropy bound", main="Huffman block coding performance", col="blue", type='l')
  
 
