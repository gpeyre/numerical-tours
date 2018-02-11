



gamma <- gamma1
displist <- round(seq(1, niter, length=10))

k <- 1

for (i in 1:niter){
  gamma <- resample( gamma + dt * normalC(gamma) )
  if (i==displist[k]){
    lw <- 1
    if (i==1 || i==niter){
      lw <- 4
    }
    if (i==1){
      cplot(gamma, type="l", pch=1, lty=1, col="red", lw=lw, add=FALSE)
    }
    else {
      cplot(gamma, type="l", pch=1, lty=1, col="red", lw=lw, add=TRUE)
    }
    k <- k+1
  }
  
  
}