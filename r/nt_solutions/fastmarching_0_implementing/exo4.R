exo4 <- function(tau,x0,x1,G){
  ####
  # Perform the full geodesic path extraction by iterating the gradient
  # descent. You must be very careful when the path become close to
  # $x_0$, because the distance function is not differentiable at this
  # point. You must stop the iteration when the path is close to $x_0$.
  ####
  
  n <- dim(G)[1]
  Geval <- function(G,x){ 
    bilinear_interpolate(G[,,1], Im(x), Re(x) ) + 1i * bilinear_interpolate(G[,,2],Im(x), Re(x)) }
  niter <- 1.5*n/tau
  # init gamma
  gamma <- c(x1)
  xtgt <- x0[1] + 1i*x0[2]
  for (i in 1:niter){
    g <- Geval(G, gamma[length(gamma)])
    gamma <- c(gamma, gamma[length(gamma)] - tau*g)
    if (abs(gamma[length(gamma)]-xtgt)<1){
      break
    }
  }
  gamma <- c(gamma, xtgt)
  return(gamma)
}



