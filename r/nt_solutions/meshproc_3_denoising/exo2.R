X1 <- X
err = c(pnoisy)

for (i in (1:12)){
  X1 <- t(tW %*% t(X1))
  err <- c(err, snr(X0,X1))
  
  if (i%%2 == 0){
    plot_mesh(data.matrix(X1), F, "grey")
  }
  if (err[length(err)] > max(err[-length(err)])){
    Xbest <- X1
  }
}