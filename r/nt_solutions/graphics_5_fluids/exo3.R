




V <- normalize(ProjI(V))
g <- f
k <- 1
niter <- 12*4

for (i in 1:niter){
  # advect
  g <- W(g,tau*U)
  V <- Wt(V,tau*U)
  # diffuse
  V <- V + tau*nu*Delta(V)
  g <- g + tau*mu*Delta(g)
  # project
  V <- ProjI(V)
  # additional constraints
  
  #display
  if (mod(i, round(niter/4)) == 0){
    imageplot(g, paste("Time =", (i*tau)), c(2,2,k))
    k <- k+1
    }
}



