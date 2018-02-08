




div <- function(Px,Py, bound="sym", order=1){
  ####
  # div - divergence operator
  # 
  # fd = div(Px,Py, options);
  # fd = div(P, options);
  # 
  # options.bound = 'per' or 'sym'
  # options.order = 1 (backward differences)
  # = 2 (centered differences)
  # 
  # Note that the -div and grad operator are adjoint
  # of each other such that 
  # <grad(f),g>=<f,-div(g)>
  #   
  #   See also: grad.
  # 
  # Copyright (c) 2007 Gabriel Peyre
  ####
  
  # retrieve number of dimensions
  nbdims <- length(dim(Px))
  
  if (nbdims >= 3){
    if (nbdims == 3){
      Py <- P[,,2]
      Px <- P[,,1]
      nbdims <- 2
    }
    else {
      Pz <- P[,,,3]
      Py <- P[,,,2]
      Px <- P[,,,1]
      nbdims <- 3
    }
  }
  
  if (bound == "sym"){
    nx <- dim(Px)[1]
    if (order == 1){
      fx <- Px - Px[c(1,1:(nx-1)),]         
      fx[1,] <- Px[1,]               # boundary
      fx[nx,] <- -Px[nx-1,] 
      
      if (nbdims>=2){
        ny <- dim(Py)[2]
        fy <- Py - Py[,c(1,(1:(ny-1)))]
        fy[,1] <- Py[,1]             # boundary
        fy[,ny] <- -Py[,ny-1]
      }
      
      if (nbdims>=3){
        nz <- dim(Pz)[3]
        fz <- Pz - Pz[,,c(1,(1:(nz-1)))]
        fz[,,1] <- Pz[,,1]           # boundary
        fy[,,nz] <- -Py[,,nz-1]
      }
    }
    
    else{
      fx <- (Px[c(2:nx,nx),] - Px[c(1,1:(nx-1)),])/2.
      fx[1,] <- + Px[2,]/2. + Px[1,] # boundary
      fx[2,] <- + Px[3,]/2. - Px[1,]
      fx[nx,] <- - Px[nx,]-Px[nx-1,]/2.
      fx[nx-1,] <- + Px[nx,]-Px[nx-2,]/2.
      
      if (nbdims>=2){
        ny <- dim(Py)[2]
        fy <- (Py[,c(2:ny,ny)] - Py[,c(1,(1:(ny-1)))])/2.
        fy[,1] <- + Py[,2]/2. + Py[,1] # boundary
        fy[,2] <- + Py[,3]/2. - Py[,1]
        fy[,ny] <- - Py[,ny]-Py[,ny-1]/2.
        fy[,ny-1] <- + Py[,ny]-Py[,ny-2]/2.
      }
      
      if (nbdims>=3){
        nz <- dim(Pz)[3]
        fz <- (Pz[,,c(2:nz,nz)] - Pz[,,c(1,(1:(nz-1)))])/2.
        fz[,,1] <- + Pz[,,2]/2. + Pz[,,1] # boundary
        fz[,,2] <- + Pz[,,3]/2. - Pz[,,1]
        fz[,,nz] <- - Pz[,,nz]-Pz[,,nz-1]/2.
        fz[,,nz-1] <- + Pz[,nz]-Pz[,nz-2]/2.
      }
    }
  }
  
  else{
    nx <- dim(Px)[1]
    if (order == 1){
      fx <- Px - Px[c(nx,1:(nx-1)),]         
      
      if (nbdims>=2){
        ny <- dim(Py)[2]
        fy <- Py - Py[,c(ny,(1:(ny-1)))]
      }
      
      if (nbdims>=3){
        nz <- dim(Pz)[3]
        fz <- Pz - Pz[,,c(nz,(1:(nz-1)))]
      }
    }
    
    else{
      nx <- dim(Px)[1]
      fx <- (Px[c(2:nx,1),] - Px[c(nx,1:(nx-1)),])
      
      if (nbdims>=2){
        ny <- dim(Py)[2]
        fy <- Py[,c(2:ny,1)] - Py[,c(1,(1:(ny-1)))]
      }
      
      if (nbdims>=3){
        nz <- dim(Pz)[3]
        fz <- Pz[,,c(2:nz,1)] - Pz[,,c(nz,(1:(nz-1)))]
      }
    }
  }
  
  # gather result
  if (nbdims == 3){
    fd <- fx+fy+fz
  }
  
  else if (nbdims == 2){
    fd <- fx+fy
  }
  
  else {
    fd <- fx
  }
  
  return(fd)
  
}


div_2 <- function(x)
{ 
  # Divergence operator
  n = dim(x)[1]
  hdiff1 = x[c(2:n, 1), , 1]
  hdiff2 = x[, c(2:n, 1), 2]
  return(hdiff1 - x[,,1] + hdiff2 - x[,,2])
}



