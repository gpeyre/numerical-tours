




circshift <- function(x, p){
  
  ####
  
  # Circular shift of an array.
  
  ##
  
  y <- as.matrix(x)
  
  if (p[1] > 1){
    
    y <- rbind( y[p[1]:dim(y)[1],,drop=F], y[1:(p[1]-1),,drop=F] )
    
    if (dim(x)[2]>0 & length(p)>1){
      
      y <- cbind( y[, p[1]:dim(y)[2],drop=F], y[, 1:(p[1]-1),drop=F] )
      
    }
    
  }
  
  else if (p[1] < -1){
    
    p[1] <- dim(y)[1] + p[1] + 2
    
    y <- rbind( y[p[1]:dim(y)[1],,drop=F], y[1:(p[1]-1),,drop=F] )
    
    if (dim(x)[2]>0 & length(p)>1){
      
      y <- cbind( y[, p[1]:dim(y)[2],drop=F], y[, 1:(p[1]-1),drop=F] )
      
    }
    
  }
  
  return(y)
  
}






circshift <- function(x, p){
  ####
  # Circular shift of an array.
  ####
  y <- as.matrix(x)
  if (p[1] > 1){
    y <- rbind( y[p[1]:dim(y)[1],,drop=F], y[1:(p[1]-1),,drop=F] )
    if (dim(x)[2]>0 & length(p)>1){
      y <- cbind( y[, p[1]:dim(y)[2],drop=F], y[, 1:(p[1]-1),drop=F] )
    }
  }
  else if (p[1] < -1){
    p[1] <- dim(y)[1] + p[1] + 2
    y <- rbind( y[p[1]:dim(y)[1],,drop=F], y[1:(p[1]-1),,drop=F] )
    if (dim(x)[2]>0 & length(p)>1){
      y <- cbind( y[, p[1]:dim(y)[2],drop=F], y[, 1:(p[1]-1),drop=F] )
    }
  }
  return(y)
}






circshift1d <- function(x, k){
  ####
  # Circularly shift a 1D vector
  ####
  return(roll(x, -k))
}


