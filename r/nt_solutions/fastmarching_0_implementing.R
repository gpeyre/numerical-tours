




exo1 <- function(x0,W){
  ####
  # Implement the Dijkstra algorithm by iterating these step while the
  # stack |I| is non empty.
  # Display from time to time the front that propagates.
  ####
  n <- dim(W)[1]
  pstart <- t(t(x0))
  DIJSTRA <- perform_dijstra_fm(W, pstart, Inf,'dijstr', 'sym',n*6)
  D <- DIJSTRA$D ; Dsvg <- DIJSTRA$Dsvg ; Ssvg <- DIJSTRA$Ssvg
  
  par(mfrow=c(2,2))
  for (i in 1:4){
    d <- Dsvg[,,i]
    d[d==Inf] <- 0

    d <- as.cimg( t(as.matrix(d)) )
    cmap_jet <- function(v){ return( rgb(v,
                                         (sin(v*2*pi)+1)/2,
                                         (cos(v*2*pi)+1)/2) ) }
    
    plot(d, colourscale=cmap_jet, interpolate = FALSE, axes = FALSE)
    
  }
  
  return(D)
  
}





exo2 <- function(x0, W){
  ####
  # Implement the FM algorithm by iterating these step while the
  # stack |I| is non empty.
  # Display from time to time the front that propagates.
  ####
  n <- dim(W)[1]
  pstart <- t(t(x0))
  DIJSTRA <- perform_dijstra_fm(W, pstart, Inf,'fm', 'sym',n*6)
  D <- DIJSTRA$D ; Dsvg <- DIJSTRA$Dsvg ; Ssvg <- DIJSTRA$Ssvg
  
  par(mfrow=c(2,2))
  for (i in 1:4){
    d <- Dsvg[,,i]
    d[d==Inf] <- 0
    
    d <- as.cimg( t(as.matrix(d)) )
    cmap_jet <- function(v){ return( rgb(v,
                                         (sin(v*2*pi)+1)/2,
                                         (cos(v*2*pi)+1)/2) ) }
    
    plot(d, colourscale=cmap_jet, interpolate = FALSE, axes = FALSE)
    
  }
  
  return(D)
  
}





exo3 <- function(x0, W){
  ####
  # Compute the distance map to these starting point using the FM algorithm.
  ####
  n <- dim(W)[1]
  pstart <- t(t(x0))
  DIJSTRA <- perform_dijstra_fm(W, pstart, Inf,'fm', 'sym',n*6)
  D <- DIJSTRA$D ; Dsvg <- DIJSTRA$Dsvg ; Ssvg <- DIJSTRA$Ssvg
  # display
  k <- 8
  displ <- function(D){ cos(2*pi*k*D/max(D) ) }
  
  cmap_jet <- function(v){ return( rgb(v,
                                       (sin(v*2*pi)+1)/2,
                                       (cos(v*2*pi)+1)/2) ) }
  
  plot(as.cimg(displ(D)), colourscale=cmap_jet, interpolate = FALSE, axes = FALSE)
  
  return(D)
  
}





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





