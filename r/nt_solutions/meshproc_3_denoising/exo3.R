Xt <- X
k <- 0
sob = c()
err = c()

for (i in (1:niter)){
  #step
  Xt <- Xt - tau * t(tL %*% t(Xt))
  #error
  err <- c(err, snr(X0,X1))
  if (i%%floor(niter/4)==0){
    k <- k+1
    plot_mesh(Xt, F)
  }
}