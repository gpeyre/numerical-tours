

my_transform <- function(x, dir){
  # my_transform - perform either FFT with energy conservation.
  # Works on array of size (w,w,a,b) on the 2 first dimensions.
  w <- dim(x)[1]
  if (dir==1) {
    y <- mvfft(x)/sqrt(w)
  } else{
    y <- mvfft(x*sqrt(w), inverse=TRUE)/length(x)
  }  
  return (y)
}


perform_stft <- function(x,w,q,n){
  #perform_stft - compute a local Fourier transform
  
  #Forward transform:
    #MF = perform_stft(M,w,q, options);
    #Backward transform:
      #M  = perform_stft(MF,w,q, options);
      
      #w is the width of the window used to perform local computation.
      #q is the spacing betwen each window.
      
      #MF(:,i) contains the spectrum around point (i-1)*q
      
      #A typical use, for an redundancy of 2 could be w=2*q+1
      
      #options.bound can be either 'per' or 'sym'
      
      #No multichannel, no Gabor stft
      
      #options.normalization can be set to
      #'tightframe': tight frame transform, with energy conservation.
      #'unit': unit norm basis vectors, usefull to do thresholding
  
  if (is.vector(x)==TRUE){
    dir = 1
  } else{
    dir = -1
  }
  
  # perform sampling
  X <- seq(from=1, to = n+1, by=q)
  
  p <- length(X)
  eta <- 1
  
  if (w%%2 == 1){
    w <- ceiling((w-1)/2)*2+1
    w1 <- floor((w-1)/2)
    dX <- seq(from=-w1, to=w1, by=1)
  } else{
    dX <- seq(from=-floor(w/2)+1,floor(w/2), by=1)
  }
  
  X1 <- matrix(rep(X,each=w),nrow=w) + t(matrix(rep((dX),each=p),nrow=p))
  #periodic boundary conditions
  X1 <- (X1-1)%%n+1 
  
  I <- X1 - 1
  
  # build a sin weight function
  W <- 0.5*(1-cos(2*pi*seq(from=0,to=w-1, by=1)/(w-1)))
  
  #renormalize the windows
  weight <- Matrix(0,nrow=1,ncol=n)
  
  for (i in (1:p)){
    weight[(I+1)[,i]] <- weight[(I+1)[,i]]+W**2
  }
  weight <- sqrt(weight)
  Weight <- t(matrix(rep((W),each=p),nrow=p))
  
  for (i in (1:p)){
    Weight[,i] <- Weight[,i]/weight[(I+1)[,i]]
  }
  
  #compute the transform
  if (dir == 1){
    y <- matrix(0, nrow=eta*w, ncol=p)
    if (w%%2== 1){
      m <- floor((eta*w+1)/2)
      w1 <- floor((w-1)/2)
      sel <- seq(from=(m-w1), to=(m+w1))-1
    } else{
       m <- floor(eta*w/2)+1
       w1 <- floor(w/2)
       sel <- seq(from=m-w1,to=m+w1-1)-1
       y[sel+1,] <- x[I+1]*Weight
    
    #perform the transform
      y <- my_transform(y,+1)
    }
  }  else{
    x <- my_transform(x,-1)  
    x <- Re(x*Weight)
    y <- Matrix(0, nrow=1, ncol=n)
    for (i in (1:p)){
      y[(I+1)[,i]] <- y[(I+1)[,i]] + x[,i]
      
    }
  }  
  return (y)  
}
