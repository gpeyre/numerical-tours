




niter <- 800
b <- array(1, c(N,N,K))
a <- b
Err_q <- rep(0, niter)

for (i in 1:niter){
  
  for (k in 1:K){
    Err_q[i] <- Err_q[i] + norm(a[,,k]*xi(b[,,k]) - P[,,k])/norm(P[,,k])
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

plot(log(Err_q), type="l", lwd = 2, col="blue", main="", xlab="", ylab="")





