


par(mfrow=c(2,2))

phi <- phi0
eps <- .Machine$double.eps
k <- 1

for (i in 1:niter){
  g0 <- grad(phi, order=2)
  d <- pmax( eps*array(1,c(n,n)), sqrt(apply(g0**2, c(1,2), sum)) )
  g <- g0/array(rep(d,2), c(dim(d),2))
  K <- d*div(g[,,1], g[,,2], order=2)
  phi <- phi + tau*K
  if (mod(i, as.integer(niter/4)) == 0){
    k <- k+1
    plot_levelset(phi, lw=2)
  }
}
