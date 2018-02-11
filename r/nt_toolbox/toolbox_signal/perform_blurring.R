




perform_blurring <- function(M, sigma, bound="sym"){
  ####
  # perform_blurring - gaussian blurs an image
  # 
  # M = perform_blurring(M, sigma, options);
  # 
  # M is the original data
  # sigma is the width of blurs (in pixels)
  # 
  # Copyright (c) 2007 Gabriel Peyre
  ####

  
  if (!(FALSE %in% (sigma == 0))){
    return(M)
  }
  
  if (length(dim(M)) > 2){
    for (i in 1:dim(M)[3]){
      M[,,i] <- perform_blurring(M[,,i], sigma, bound)
    }
  }
  
  n <- max(dim(M))
  
  eta <- 4
  p <- round((sigma*eta)/2.)*2+1
  p <- pmin(p,(round(n/2.)*2-1)*rep(1,length=length(p)))
  
  A <- c(1., 1.)
  if (length(dim(M))==0){
    A <- 1 #1D
  }
  
  h <- compute_gaussian_filter(p*A,sigma/(4.*n),n*A)
  M <- perform_convolution(M, h, bound)
  
  return(M)

}





compute_gaussian_filter <- function(n,s,N){
  ####
  # compute_gaussian_filter - compute a 1D or 2D Gaussian filter.
  # 
  # f = compute_gaussian_filter(n,s,N);
  # 
  # 'n' is the size of the filter, odd for no phase in the filter.
  # (if too small it will alterate the filter).
  # use n=[n1,n2] for a 2D filter or n = [n1] for a 1D filter
  # 's' is the standard deviation of the filter.
  # 'N' is the size of the big signal/image (supposed to lie in [0,1] or [0,1]x[0,1]).
  # use N=[N1,N2] for a 2D filter or N = [N1] for a 1D filter
  # 
  # The equation (in 1D) is
  # f[k] = exp( -(x(k)^2/(2*s^2)) );
  # where x spans [-1/2,1/2].
  # 
  # The filter is normalised so that it sums to 1.
  # 
  # Copyright (c) 2004 Gabriel Peyre
  ####
  
  nd <- 1
  if (length(n)>1 & n[2]>1){
    nd <- 2
  }
  
  if (nd==2 & length(s)==1){
    s <- c(s,s)
  }
  
  if (nd==2 & length(N)==1){
    N <- c(N,N)
  }
  
  if (nd==1){
    f <- build_gaussian_filter_1d(n,s,N)
  }
  
  else {
    f <- build_gaussian_filter_2d(n,s,N)
  }
  
  return(f)

}





build_gaussian_filter_2d <- function(n,s,N=c()){
  ####
  # build_gaussian_filter_2d - compute a 2D Gaussian filter.
  # 
  # f = build_gaussian_filter_2d(n,s,N);
  # 
  # 'n' is the size of the filter, odd for no phase in the filter.
  # (if too small it will alterate the filter).
  # 's' is the standard deviation of the filter.
  # 'N' is the size of the big image (supposed to lie in [0,1]x[0,1]).
  # 
  # The filter is normalised so that it sums to 1.
  # 
  # Copyright (c) 2004 Gabriel Peyre
  ####
  
  n <- as.vector(n)
  s <- as.vector(s)
  N <- as.vector(N)
  
  if (length(N)==0){
    N <- n
  }
  
  if (length(N)==1 | N[1]==1){
    N <- c(N, N)
  }
  
  if (length(s)==1 | s[1]==1){
    s <- c(s, s)
  }
  
  if (length(s[s<=0]) > 0){
    f <- array(0, n)
    f[round((n-1)/2)] <- 1 # +1 ?
    return(f)
  }
  
  x <- (0:(n[1]-1) - (n[1]-1)/2.)/(N[1]-1)
  y <- (0:(n[2]-1) - (n[2]-1)/2.)/(N[2]-1)
  grid <- meshgrid_2d(y, x)
  Y <- grid$X ; X <- grid$Y
  f <- exp( -(X**2/(2*s[1]**2)) - (Y**2/(2*s[2]**2)) )
  f <- f / sum(f)
  return(f)
}




build_gaussian_filter_1d <-function(n,s,N=c()){
  ####
  # build_gaussian_filter_1d - compute a Gaussian filter.
  # 
  # f = build_gaussian_filter_1d(n,s,N);
  # 
  # Copyright (c) 2004 Gabriel Peyre
  ####
  
  if (length(N)==0){
    N <- n
  }
  
  n <- n[1]
  s <- s[1]
  N <- N[1]
  
  if (s <= 0){
    f <- rep(0, n)
    f[round((n-1)/2)] <- 1
    return(f)
  }
  
  x <- ( 0:(n-1) - (n-1)/2.)/(N-1)
  f <- exp( -x**2 / (2*s**2) )
  f <- f / sum(f)
  return(f)
}






perform_convolution <- function(x,h,bound="sym"){
  ####
  # perform_convolution - compute convolution with centered filter.
  # 
  # y = perform_convolution(x,h,bound);
  # 
  # The filter 'h' is centred at 0 for odd
  # length of the filter, and at 1/2 otherwise.
  # 
  # This works either for 1D or 2D convolution.
  # For 2D the matrix have to be square.
  # 
  # 'bound' is either 'per' (periodic extension) 
  # or 'sym' (symmetric extension).
  # 
  # Copyright (c) 2004 Gabriel Peyre
  ####
  
  if (!(bound %in% c("sym", "per"))){
    warning( 'bound should be sym or per' )
    return()
  }
  
  if (length(dim(x))==3 & dim(x)[3]<4){
    #for color images
    y <- x
    for (i in 1:dim(x)[3]){
      y[,,i] <- perform_convolution(x[,,i],h, bound)
    }
    return(y)
  }
  
  if (length(dim(x))==3 & dim(x)[3]>=4){
    warning( 'Not yet implemented for 3D array, use smooth3 instead.' )
    return()
  }
  
  n <- dim(x)
  if (length(n)==0){ n <- length(x)}
  p <- dim(h)
  if (length(p)==0){ p <- length(h)}
  
  nd <- length(dim(x))
  if (nd == 0){ nd <- 1}
  
  if (nd == 1){
    n <- length(x)
    p <- length(h)
  }
  
  if (bound == "sym"){
    
    #################################
    # symmetric boundary conditions #
    d1 <- as.integer(p)/2   # padding before
    d2 <- p - d1 - 1        # padding after
    
    if (nd==1){
      ############# 1D #############
      nx <- length(x)
      xx <- array(c(x[(length(x)-1):(d1+1)],
                    x,
                    x[(nx-d2-1):nx]),
                  c(3, length(x)))
      y <- convolve(xx, h)##############
      y <- y[p+1:nx-p-1]
    }
    else if (nd==2){
      ############# 2D #############
      #double symmetry
      nx <- dim(x)[1] ; ny <- dim(x)[2]
      xx <- x
      
      top <- xx[(nx-1):(d1[1]+1),]
      bottom <- xx[(nx-d2[1]-1):nx,]
      temp <- array(0, c(dim(top)[1]+dim(xx)[1]+dim(bottom)[1], ny))
      temp[1:dim(top)[1],] <- top
      temp[(dim(top)[1]+1):(dim(top)[1]+1+dim(xx[1])),] <- xx
      temp[(dim(top)[1]+1+dim(xx[1])):dim(temp)[1],] <- bottom
      xx <- temp
      
      left <- xx[,(ny-1):(d1[2]+1)]
      right <- xx[,(ny-d2[2]-1):ny]
      temp <- array(0, c(dim(xx)[1], dim(left)[2]+dim(xx)[2]+dim(right)[2]))
      temp[,(1:dim(left)[2])] <- left
      temp[,(1+dim(left)[2]):(1+dim(left)[2]+dim(xx)[2])] <- xx
      temp[,(1+dim(left)[2]+dim(xx)[2]):dim(temp)[2]] <- right
      
      y <- convolve(xx, h)
      y <- y[(2*d1[1]+1):(2*d1[1]+n[1]+1), (2*d1[2]+1):(2*d1[2]+n[2]+1)]  
    }
  }
  
  else{
    
    ################################
    # periodic boundary conditions #
    
    if (p > n){
      print(paste("p =",p))
      print(paste("n =",n))
      warning('h filter should be shorter than x.')
      return()
    }
    print(p)
    n <- as.vector(n)
    p <- as.vector(p)
    d <- as.integer((p-1)/2)
    if (nd == 1){
      h <- c(h[(d+1):length(h)],
             rep(0, n-p),
             h[1:d])
      fft_mult <- fft(x)*fft(h)
      y <- Re(fft(fft_mult, inverse=TRUE)/length(fft_mult))
    }
    else{
      top <- h[(as.integer((d[1]))+1):dim(h)[1],]
      middle <- array(0, c(n[1]-p[1], p[2]))
      bottom <- h[(1:as.integer((d[1]))),]
      temp <- array(0, c(dim(top)[1]+dim(middle)[1]+dim(bottom)[1], dim(h)[2]))
      temp[1:dim(top)[1],] <- top
      temp[(dim(top)[1]+1):(dim(top)[1]+dim(middle)[1]),] <- middle
      temp[(1+dim(top)[1]+dim(middle)[1]):dim(temp)[1],] <- bottom
      h <- temp
      
      
      left <- h[,(as.integer((d[2]))+1):dim(h)[2]]
      middle <- array(0, c(n[1], n[2]-p[2]))
      right <- h[,(1:as.integer((d[2])))]
      temp <- array(0, c(dim(h)[1], dim(left)[2]+dim(middle)[2]+dim(right)[2]))
      temp[,1:dim(left)[2]] <- left
      temp[,(dim(left)[2]+1):(dim(left)[2]+dim(middle)[2])] <- middle
      temp[,(1+dim(left)[2]+dim(middle)[2]):dim(temp)[2]] <- right
      h <- temp
      
      fft_mult <- fft(x)*fft(h)
      y <- Re(fft(fft_mult, inverse=TRUE)/length(fft_mult))
    }
    return(y)
    
    
    
  }
  
  
  
}








