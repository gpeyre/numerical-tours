



slist <- c(4, 6, 10, 15)

for (i in 1:length(slist)){
  sigma <- slist[i]
  
  g <- grad(blur(f0, sigma))
  H <- hessian(blur(f0, sigma))
  a <- H[,,c(1,2)] * array(rep(g[,,1],2), c(n,n,2)) + H[,,c(2,3)] * array(rep(g[,,2],2), c(n,n,2))
  plot_levelset(apply(a*g, c(1,2), sum), f0, paste("sigma =", sigma), sbpt=c(2,2,i))
  
}