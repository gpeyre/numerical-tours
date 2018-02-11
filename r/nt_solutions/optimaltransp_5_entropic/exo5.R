



# Insert your code here.


m <- 5 
grid <- meshgrid_2d(seq(0,1,length=m), seq(0,1,length=m))
T <- grid$X ; S <- grid$Y
T <- as.vector(T)
S <- as.vector(S)
niter <- 25 #1000

for (j in 1:(m**2)){
  # weights
  lambd <- c(S[j]*T[j], (1-S[j])*T[j], S[j]*(1-T[j]), (1-S[j])*(1-T[j]))
  # computation
  b <- array(1, c(N,N,K))
  a <- b
  
  for (i in 1:niter){
    
    for (k in 1:K){
      a[,,k] <- P[,,k]/xi(b[,,k])
    }
    
    q <- rep(0, N)
    
    for (k in 1:K){
      q <- q + lambd[k] * log(pmax(1e-19*rep(1, length(b[,,k])), b[,,k]*xi(a[,,k])))
    }
    
    q <- exp(q)
    
    for (k in 1:K){
      b[,,k] <- q/xi(a[,,k])
    }
    
  }
  
  imageplot(as.cimg(q), '', c(m,m,j))
  
}