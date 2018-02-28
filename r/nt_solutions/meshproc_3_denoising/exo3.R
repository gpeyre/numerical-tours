Xt <- X
k <- 0
sob = c()
err = c()

for (i in (1:niter)){
  #step
  Xt <- Xt - tau * t(tL %*% t(Xt))
  #error
  err <- c(err, snr(X0,Xt))
  if (i%%floor(niter/4)==0){
    k <- k+1
    trimesh(t(F+1),data.matrix(t(Xt)),main=paste('T=',Tmax*k/4), col="grey")
  }
}