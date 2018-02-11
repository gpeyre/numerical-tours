



perform_wavelet_transf <- function(f, Jmin, dir, filter = "9-7",separable = 0, ti = 0){
  ####
  # perform_wavelet_transf - peform fast lifting transform
  # 
  # y = perform_wavelet_transf(x, Jmin, dir, filter = "9-7",separable = 0, ti = 0);
  # 
  # Implement 1D and 2D symmetric wavelets with symmetric boundary treatements, using
  # a lifting implementation.
  # 
  # filter gives the coefficients of the lifting filter.
  # You can use h='linear' or h='7-9' to select automatically biorthogonal 
  # transform with 2 and 4 vanishing moments.
  # 
  # You can set ti=1 to compute a translation invariant wavelet transform.
  # 
  # You can set separable=1 to compute a separable 2D wavelet
  # transform.
  # 
  # Copyright (c) 2008 Gabriel Peyre
  ####
  
  #copy f
  x <- f

  #convert Jmin to int
  Jmin = as.integer(Jmin)
  
  # detect dimensionality 
  d <- length(dim(x))
  if (d==4){
    x <- as.matrix(x) # 'cimg' image
    d <- 2}
  if (max(d,0)==0){d <- 1} # vector type
  # P/U/P/U/etc the last coefficient is scaling
  if (filter %in% c("linear","5-3")){
    h <- c(1/2, 1/4, np.sqrt(2))
  }
  
  else if (filter %in% c("9-7","7-9")){
    h <- c(1.586134342, -0.05298011854, -0.8829110762, 0.4435068522, 1.149604398)
  }
  
  else {
    warning("Unknown filter")
  }
  
  if ( d==2 & separable==1 ){
    ti <- 0
    if (ti==1){
      warning("Separable does not works for translation invariant transform")
    }
    
    # perform a separable wavelet transform
    n <- dim(x)[1]

    if (dir==1){
      for (i in 1:n){
        x[,i] <- perform_wavelet_transf(x[,i], Jmin, dir, filter, separable, ti)
      }
      for (i in 1:n){
        x[i,] <- t(perform_wavelet_transf(t(x[i,]), Jmin, dir, filter, separable, ti))
      }
    }
    else{
      for (i in 1:n){
        x[i,] <- t(perform_wavelet_transf(t(x[i,]), Jmin, dir, filter, separable, ti))
      }
      for (i in 1:n){
        x[,i] <- perform_wavelet_transf(x[,i], Jmin, dir, filter, separable, ti)
      }
    }
  }
  
  # number of lifting steps
  if (d==1){
    n <- length(x)
  }
  else{
    n <- dim(x)[2]
  }
  m <- (length(h)-1)%/%2
  

  Jmax <- as.integer(log2(n)-1)
  jlist <- Jmax:Jmin
  
  if (dir==-1){
    jlist <- Jmin:Jmax
  }
  
  if (ti==0){
    # subsampled
    for (j in jlist){
      if (d==1){
        x[1:(2**(j+1))] <- lifting_step(x[1:2**(j+1), drop=F], h, dir)
      }
      else{
        x[1:(2**(j+1)),1:(2**(j+1))] <- lifting_step(x[1:(2**(j+1)),1:(2**(j+1))], h, dir) 
        x[1:(2**(j+1)),1:(2**(j+1))] <- t(lifting_step(t(x[1:(2**(j+1)),1:(2**(j+1))]), h, dir)) 
        
      }
    }
  }
  
  else{
    # TI
    nJ <- Jmax - Jmin + 1
    if (dir==1 & d==1){
      x <- array(rep(x, each=(nJ + 1)), c((nJ + 1), 1, length(x)))
    }
    else if (dir==1 & d==2){
      x <- array(rep(x, each=(3*nJ+1)), c((3*nJ + 1), dim(x)[1], dim(x)[2]))
    }
    
    for (j in jlist){
      dist <- 2**(Jmax - j)
      
      if (d==1){
        if (dir==1){
          x[1:(j-Jmin+2)] <- lifting_step_ti(x[1], h, dir, dist)
          }
        else{
          x[1] <- lifting_step_ti(x[1:(j-Jmin+2)], h, dir, dist)
          }
      }
      
      else{
        dj <- 3*(j-Jmin)
        
        if (dir==1){
          if (length(dim(x))==2){
            x[c(1,(dj+2)),] <- lifting_step_ti(x[1,], h, dir, dist)

            x[c(1,(dj+3)),] <- lifting_step_ti(t(x[1,]), h, dir, dist)
            x[1,] <- t(x[1,])
            x[(dj+3),] <- t(x[(dj+3),])
            
            x[c((2+dj),(4+dj)),] <- lifting_step_ti(t(x[(dj+2),]), h, dir, dist)
            x[(dj+2),] <- t(x[(dj+2),])
            x[(dj+4),] <- t(x[(dj+4),])
            
            }
          else if (length(dim(x))==3){
            x[c(1,(dj+2)),,] <- lifting_step_ti(x[1,,], h, dir, dist)
            
            x[c(1,(dj+3)),,] <- lifting_step_ti(t(x[1,,]), h, dir, dist)
            x[1,,] <- t(x[1,,]) 
            x[(dj+3),,] <- t(x[(dj+3),,]) 
            
            x[c((2+dj),(4+dj)),,] <- lifting_step_ti(t(x[(dj+2),,]), h, dir, dist)
            x[(dj+2),,] <- t(x[(dj+2),,])
            x[(dj+4),,] <- t(x[(dj+4),,])
            
            }
      
        }
        
        else{
     
          x[dj+2,,] <- t(x[dj+2,,])
          x[dj+4,,] <- t(x[dj+4,,])
          
          x[dj+2,,] <- t((lifting_step_ti(x[c(2+dj, 4+dj),,], h, dir, dist)))
          
          x[1,,] <- t(x[1,,])
          x[dj+3,,] <- t(x[dj+3,,])
          x[1,,] <- t((lifting_step_ti(x[c(1,dj+3),,], h, dir, dist)))
          
          x[1,,] <- lifting_step_ti(x[c(1,dj+2),,], h, dir, dist) 
          
        }
        
      }
   
    }
  
    if (dir==-1){
      if (length(dim(x))==3){x <- x[1,,]}
      else if (length(dim(x))==2){x <- x[1,]}
      else if (length(dim(x))==1){x <- x[1]}
    }
    
  }
  
  return(x)
  
}


###########################################################################
###########################################################################
###########################################################################



lifting_step <- function(x0, h, dir){
  
  #copy x
  x <- x0
  
  # number of lifting steps
  m <- (length(h) - 1)%/%2
  
  if (dir==1){
    # split
    d <- x[seq(2,dim(x)[1],2),, drop=F]
    x <- x[seq(1,dim(x)[1],2),, drop=F]
    for (i in (0:(m-1))){
      d <- d - h[2*i+1] * (x + rbind(x[2:dim(x)[1],, drop=F], x[dim(x)[1],, drop=F]))
      x <- x + h[2*i+2] * (d + rbind(d[1,, drop=F],d[1:(dim(d)[1]-1),, drop=F]))
    }
    x <- rbind(x*h[length(h)], d/h[length(h)])
  }
  
  else{
    # retrieve detail coefs
    end <- dim(x)[1]
    if (max(end,0)==0){end <- length(x)}
    d <- x[(as.integer(end/2) + 1):dim(x)[1],, drop=F]*h[length(h)]
    x <- x[1:as.integer(end/2),, drop=F]/h[length(h)]
    for (i in (m:1)){
      x <- x - h[2*i] * (d + rbind(d[1,, drop=F],d[1:(dim(d)[1]-1),, drop=F]))
      d <- d + h[2*i-1] * (x + rbind(x[2:dim(x)[1],, drop=F], x[dim(x)[1],, drop=F]))
    }
    # merge
    x1 <- rbind(x,x)
    x1[seq(1,dim(x1)[1],2),] <- x
    x1[seq(2,dim(x1)[1],2),] <- d
    x <- x1
  }
  
  return(x)
  
}


###########################################################################
###########################################################################
###########################################################################




lifting_step_ti <- function(x0, h, dir, dist){
  
  #copy x
  #if (length(dim(x0))==2){
  #  x <- x0[,, drop=F]
  #}
  #else if (length(dim(x0))==3){
  #  x <- x0[,,, drop=F]
  #}
  x <- x0

  # number of lifting steps
  m <- (length(h) - 1)%/%2
  if (length(dim(x))==2){
    n <- dim(x[1,,drop=F])[2]
  }
  else if (length(dim(x))==3){
    n <- dim(x[1,,,drop=F])[2]
    }
  else{
    warning("Not implemented")
  }
  
  
  s1 <- 1:n + dist
  s2 <- 1:n - dist
  
  # boundary conditions
  s1[s1 > n] <- 2*n - s1[s1 > n]
  s1[s1 < 1] <- 2   - s1[s1 < 1]
  
  s2[s2 > n] <- 2*n - s2[s2 > n]
  s2[s2 < 1] <- 2   - s2[s2 < 1]
  
  if (dir == 1){
    # split
    if (length(dim(x)) == 2){
      d <- array(x, c(1, dim(x)[1], dim(x)[2]))
    }
    else{ d <- x }
    for (i in 0:(m-1)){
      if (length(dim(x)) == 2){
        x <- array(x, c(1, dim(x)[1], dim(x)[2]))
      }

      d <- d - h[2*i+1] * (x[,s1,,drop=F ] + x[,s2,,drop=F ])
      x <- x + h[2*i+2] * (d[,s1,,drop=F ] + d[,s2,,drop=F ])

    }

    #merge
    if (length(dim(x))==2){ x <- rbind((x*h[length(h)]), (d/h[length(h)]))}
    else if (length(dim(x))==3){
      temp <- array( rep(0, length(x[1,,])+length(d[1,,])), c(dim(x)[1]+dim(d)[1], dim(x)[2], dim(x)[3]))
      temp[1:dim(x)[1],,] <- x*h[length(h)]
      temp[(dim(x)[1]+1):(dim(x)[1]+dim(d)[1]),,] <- d/h[length(h)]
      x <- temp
     }
  }
  
  else{
    # retrieve detail coefs
    if (length(dim(x))==2){
      d <- x[2,]*h[length(h)]
      x <- x[1,]/h[length(h)]
    }
    else if (length(dim(x))==3){
      d <- x[2,,]*h[length(h)]
      x <- x[1,,]/h[length(h)]
    }
    
    
    for (i in m:1){
      x <- x - h[2*i]*(d[s1,] + d[s2,]) 
      d <- d + h[2*i-1]*(x[s1,] + x[s2,])
      
    }
    # merge
    x <- (x + d)/2
  }
  
  return(x)
  
}







