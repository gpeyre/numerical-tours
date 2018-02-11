



grad <- function(M, bound="sym", order=1){
  ####
  # grad - gradient, forward differences
  # 
  # [gx,gy] = grad(M, options);
  # or
  # g = grad(M, options);
  # 
  # options.bound = 'per' or 'sym'
  # options.order = 1 (backward differences)
  # = 2 (centered differences)
  # 
  # Works also for 3D array.
  # Assme that the function is evenly sampled with sampling step 1.
  # 
  # See also: div.
  # 
  # Copyright (c) Gabriel Peyre
  ####
  
  
  
  # retrieve number of dimensions
  nbdims <- length(dim(M))
  if (nbdims==0){nbdims <- 1}
  
  if (bound == "sym"){
    nx <- dim(M)[1]
    if (order == 1){
      fx <- M[c(2:nx, nx),] - M
    }
    else{
      fx <- (M[c(2:nx, nx),] - M[c(1, 1:(nx-1)),])/2
      # boundary
      fx[1,] <- M[2,]-M[1,]
      fx[nx,] <- M[nx,]-M[nx-1,]
    }
    
    if (nbdims >= 2){
      ny <- dim(M)[2]
      if (order == 1){
        fy <- M[,c(2:ny, ny)] - M
      }
      else{
        fy <- (M[,c(2:ny, ny)] - M[,c(1,1:(ny-1))])/2
        # boundary
        fy[,1] <- M[,2]-M[,1]
        fy[,ny] <- M[,ny]-M[,ny-1]
      }
    }
    
    if (nbdims >= 3){
      nz <- dim(M)[3]
      if (order == 1){
        fz <- M[,,c(2:nz, nz)] - M
      }
      else{
        fz <- (M[,c(2:nz, nz)] - M[,c(1,1:(nz-1))])/2
        # boundary
        fz[,,1] <- M[,,2]-M[,,1]
        fz[,,nz] <- M[,,nz]-M[,,nz-1]
      }
    }
  }
  
  else{
    nx <- dim(M)[1]
    if (order == 1){
      fx <- M[c(2:nx, 1),] - M
    }
    else{
      fx <- (M[c(2:nx, 1),] - M[c(nx, 1:(nx-1)),])/2
    }
    
    if (nbdims >= 2){
      ny <- dim(M)[2]
      if (order == 1){
        fy <- M[,c(2:ny, 1)] - M
      }
      else{
        fy <- (M[,c(2:ny, 1)] - M[,c(ny,1:(ny-1))])/2
      }
    }
    
    if (nbdims >= 3){
      nz <- dim(M)[3]
      if (order == 1){
        fz <- M[,,c(2:nz, 1)] - M
      }
      else{
        fz <- (M[,c(2:nz, 1)] - M[,c(nz,1:(nz-1))])/2
      }
    }
  }
  
  if (nbdims == 2){
    out <- array(rep(0, length=length(M)*2), c(dim(M),2))
    out[,,1] <- fx ; out[,,2] <- fy
  }
  
  else if (nbdims == 3){
    out <- array(rep(0, length=length(M)*3), c(dim(M),3))
    out[,,,1] <- fx ; out[,,,2] <- fy ; out[,,,3] <- fz
  }
  
  return(out)

}


grad_2 <- function(x){
  n = dim(x)[1]
  hdiff = x[c((2:n),1),] - x[,]
  vdiff = x[,c((2:n),1)] - x[,]
  return (array(c(hdiff, vdiff), dim=c(n, n, 2)))
}

grad_3 = function (x) 
{
    n = dim(x)[1]
    hdiff = x[, ] - x[c(n, 1:n - 1), ]
    vdiff = x[, ] - x[, c(n, 1:n - 1)]
    return(array(c(hdiff, vdiff), dim = c(n, n, 2)))
}


