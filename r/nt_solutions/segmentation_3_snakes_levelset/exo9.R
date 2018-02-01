




par(mfrow=c(2,2))

phi <- phi0
k <- 1

for (i in 1:niter){
  gD <- grad(phi, order=2)
  d <- pmax(eps*array(1, c(n,n)), sqrt(apply(gD**2, c(1,2), sum)))
  g <- gD/array(rep(d,2), c(dim(d),2))
  G <- d*div(g[,,1], g[,,2], order=2) - lambd*(f0-c1)**2 + lambd*(f0-c2)**2
  phi <- phi + tau*G
  if (mod(i, as.integer(niter/4))==0){
    k <- k+1
    plot_levelset(phi, f0, lw=2)
  }
  
  
}