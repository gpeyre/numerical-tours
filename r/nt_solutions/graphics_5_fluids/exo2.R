




rho <- .25
niter <- 12*4
k <- 1
f1 <- f

for (i in 2:niter){
  f1 <- W(f1, rho*U)
  if (mod(i, round(niter/4)) == 0){
    imageplot(f1, paste("t =", (i*rho)), c(2,2,k))
    k <- k+1
  }
}
