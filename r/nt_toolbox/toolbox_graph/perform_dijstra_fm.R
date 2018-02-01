




perform_dijstra_fm <- function(W, pstart, niter=Inf, method='dijstr', bound='sym', svg_rate=10){
  
  ####
  #   perform_fm_dijstra - *slow* (matlab) implementation of Dijstra and FM
  #
  #   [D,Dsvg,Ssvg] = perform_fm_dijstra(W, pstart, options);
  #
  #   W is an (n,n) metric matrix.
  #   pstart is a (2,k) starting points.
  #   options.method is either 'fm' or 'dijstra'
  #   
  #   D is the final distance map to pstart
  #   options.svg_rate gives the rate at wich Dsvg and Ssvg is filled.
  #   options.niter can be used to limit the total number of steps (partial propagation). 
  #   
  #   Copyright (c) 2014 Gabriel Peyre
  ####
  
  ##
  # Size.
  
  n <- dim(W)[1]
  
  
  ##
  # The four displacement vector to go to the four neightbors.
  
  neigh <- array(c(1, 0, -1, 0, 0, 1, 0, -1), c(2,4))
  
  
  ##
  # For simplicity of implementation, we use periodic boundary conditions.
  
  # boundary = @(x)x.*(x<=n & x>0) + (2-x).*(x<=0) + (2*n-x).*(x>n);
  
  symmetrize <- function(x,n){
    if (x<=0){
      x <- 2-x
    }
    else if (x>n){
      x <- 2*n-x
    }
    return(x)
  }
  
  
  if (bound=='per'){
    boundary <- function(x){ mod(x,n) }
  }
  else{
    boundary <- function(x){ c(symmetrize(x[1],n), symmetrize(x[2],n)) } # todo
  }
  
  
  ##
  # For a given grid index |k|, and a given neighboring index k in \({1,2,3,4}\), 
  # |Neigh(k,i)| gives the corresponding grid neigboring index.
  
  ind2sub1 <- function(k){ c(as.integer( (k-1-mod(k-1,n))/n + 1 ), mod(k-1,n) + 1) }
  sub2ind1 <- function(u){ as.integer( (u[1]-1)*n + (u[2]-1) + 1 ) }
  Neigh <- function(k,i){ sub2ind1(boundary(ind2sub1(k) + neigh[,i])) }
  extract   <- function(x,I){ x[I] }
  extract1d <- function(x,I){ extract(as.vector(t(x)),I) }
  
  ##
  # Stack of starting points.
  
  nstart <- dim(pstart)[2]
  I <- rep(0, nstart)
  for (i in 1:nstart){
    I[i] <- as.integer( sub2ind1(pstart[,i]) )
  }
  
  
  ##
  # Initialize the distance to \(+\infty\), excepted for the boundary conditions.
  
  D <- array(0, c(n,n)) + Inf # current distance
  for (i in 1:nstart){
    D[pstart[1,i],pstart[2,i]] <- 0
  }
  
  
  ##
  # Initialize the state to 0 (unexplored), excepted for the boundary point to \(1\)
  # (front).
  
  S <- array(0, c(n,n))
  for (i in 1:nstart){
    S[pstart[1,i],pstart[2,i]] <- 1 # open
  }
  
  ##
  # Run!
  
  iter <- 0
  q = 100  # maximum number of saves
  Dsvg <- array( 0, c(n,n,q) )
  Ssvg <- array( 0, c(n,n,q) )
  while ( length(I)!=0 & (iter<=niter) ){
    # print(length(I)!=0 & (iter<=niter))      
    iter <- iter+1
    # print(length(I))
    if (iter==niter){
      break
    }
    # pop from stack
    j <- order( extract1d(D,I)  )
    j <- j[1]
    i <- I[j]       
    a <-  I[j]
    I <- I[-c(j)]
    # declare dead
    u <- ind2sub1(i)
    S[u[1],u[2]] <- -1
    # Make a list of neighbors that are not dead
    J <- c()
    for (k in 1:4){
      j <- Neigh(i,k)
      if ( extract1d(S,j)!=-1 ){
        # add to the list of point to update
        J <- c(J, j) 
        if ( extract1d(S,j)==0 ){
          # add to the front
          u <- ind2sub1(j)
          S[u[1],u[2]] <- 1
          I <- c(I, j)
        }
      }
    }
    # update neighbor values
    DNeigh <- function(D,k){ extract1d(D,Neigh(j,k)) }
    for (j in J){
      dx <- min(DNeigh(D,1), DNeigh(D,2))
      dy <- min(DNeigh(D,3), DNeigh(D,4))
      u <- ind2sub1(j)
      w <- extract1d(W,j);
      if (method == "dijstr"){
        D[u[1],u[2]] <- min(dx + w, dy + w)
      }
      else{
        Delta <- 2*w - (dx-dy)**2
        if (Delta>=0){
          D[u[1],u[2]] <- (dx + dy + sqrt(Delta))/ 2
        }
        else{
          D[u[1],u[2]] <- min(dx + w, dy + w)
        }
      }
    }
    # svd
    t <- iter/svg_rate
    if ((mod(iter,svg_rate)==0) & (t<q)){
      Dsvg[,,t] <- D
      Ssvg[,,t] <- S
    }
    
  }
  
  Dsvg <- Dsvg[,,1:t]
  Ssvg = Ssvg[,,1:t]
  return(list(D = D, Dsvg = Dsvg, Ssvg = Ssvg))
  
}