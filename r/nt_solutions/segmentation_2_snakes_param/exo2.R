


gamma <- gamma0
displist <- round(seq(1, niter, length=10))
k <- 1
imageplot(t(W))
for (i in 1:niter){
  N <- normal(gamma)
  g <- EvalW(gamma)*normalC(gamma) - N*dotp(EvalG(gamma), N)
  gamma <- resample( gamma + dt*g )
  if (i==displist[k]){
    lw <- 1
    if (i==1 || i==niter){
      lw <- 4
    }
    cplot(gamma, type="l", pch=1, lw=lw, lty=1, col="red", add=TRUE)
    k <- k+1
  }
}