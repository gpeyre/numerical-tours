



grid <- meshgrid_2d(1:n, 1:n)
Y <- grid$X ; X <- grid$Y
r <- n/3
c <- c(n,n)/2
phi0 <- pmax(abs(X-c[1]), abs(Y-c[2])) - r