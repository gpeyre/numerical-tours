



# Functions


norm <- function(v){
  ####
  # Euclidian norm of a vector
  ####
  return(sqrt(sum(v**2)))
}



roll <- function(x, n, axis=1){
  ####
  # Roll of a vector or matrix (column roll in case of matrix)
  ####
  if(n==0){
    return(x)
  }
  else{
    if (length(dim(x)) <= 1){
      return(c( tail(x,n) , head(x,-n) ))
    }
    else if (length(dim(x)) == 2){
      return( apply(x,2,roll,n) )
    }
    else if (length(dim(x)) == 3){
      if (axis==1){
        if (n<0){ n <- (dim(x)[1]+n) }
        a <- x[1:(dim(x)[1]-n),,, drop=F]
        b <- x[(dim(x)[1]-n+1):(dim(x)[1]),,, drop=F]
        c <- array(0, dim(x))
        c[1:n,,] <- b
        c[(n+1):dim(x)[1],,] <- a
        return(c)
      }
      if (axis==2){
        if (n<0){ n <- dim(x)[2]+n }
        a <- x[,1:(dim(x)[2]-n),, drop=F]
        b <- x[,(dim(x)[2]-n+1):(dim(x)[2]),, drop=F]
        c <- array(0, dim(x))
        c[,1:n,] <- b
        c[,(n+1):dim(x)[2],] <- a
        return(c)
      }
      if (axis==3){
        if (n<0){ n <- dim(x)[3]+n }
        a <- x[,,1:(dim(x)[3]-n), drop=F]
        b <- x[,,(dim(x)[3]-n+1):(dim(x)[3]), drop=F]
        c <- array(0, dim(x))
        c[,,1:n] <- b
        c[,,(n+1):dim(x)[3]] <- a
        return(c)
      }
    }
    else{
      warning('Not implemented')
    }
  }
  
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



reverse <- function(x){
  ####
  # Reverse a vector
  ####
  return(rev(x))
}



rescale <- function(f, a=0, b=1){
  ####
  # Rescale linearly the dynamic of a vector to fit within a range [a,b]
  ####
  v <- max(f) - min(f)
  g <- (f - min(f))
  if (v>0){
    g <- g/v
  }
  return(a + g*(b-a))
}

clamp = function(x, a=0, b=1)
{
    
    "clamp - clamp a value

       y = clamp(x,a,b);

     Default is [a,b]=[0,1].
    "
    return (pmin(pmax(x,a),b))
}





ravel <- function(M){
  ####
  # return the 1D-array correponding to matrix M (row-wise)
  ####
  return( as.vector(t(M)) )
}





