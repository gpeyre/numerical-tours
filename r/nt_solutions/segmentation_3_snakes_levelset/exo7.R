




grid <- meshgrid_2d(1:n, 1:n)
Y <- grid$X ; X <- grid$Y
k <- 4 #number of circles
r <- .3*n/k
phi0 <- array(Inf, c(n,n))

for (i in 1:k){
  for (j in 1:k){
    c <- (c(i,j) - 1)*(n/k) + (n/k)*.5
    phi0 <- pmin(phi0,
                 sqrt(abs(X-c[1])**2 + abs(Y-c[2])**2) - r)
  }
}

par(mfrow=c(1,2))

plot_levelset(phi0,lw=2)

plot_levelset(phi0, f0, lw=2)