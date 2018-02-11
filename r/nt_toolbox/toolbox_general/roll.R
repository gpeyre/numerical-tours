





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
    
    else if (length(dim(x)) == 4){
      if (axis==1){
        if (n<0){ n <- (dim(x)[1]+n) }
        a <- x[1:(dim(x)[1]-n),,,, drop=F]
        b <- x[(dim(x)[1]-n+1):(dim(x)[1]),,,, drop=F]
        c <- array(0, dim(x))
        c[1:n,,,] <- b
        c[(n+1):dim(x)[1],,,] <- a
        return(c)
      }
      if (axis==2){
        if (n<0){ n <- dim(x)[2]+n }
        a <- x[,1:(dim(x)[2]-n),,, drop=F]
        b <- x[,(dim(x)[2]-n+1):(dim(x)[2]),,, drop=F]
        c <- array(0, dim(x))
        c[,1:n,,] <- b
        c[,(n+1):dim(x)[2],,] <- a
        return(c)
      }
      if (axis==3){
        if (n<0){ n <- dim(x)[3]+n }
        a <- x[,,1:(dim(x)[3]-n),, drop=F]
        b <- x[,,(dim(x)[3]-n+1):(dim(x)[3]),, drop=F]
        c <- array(0, dim(x))
        c[,,1:n,] <- b
        c[,,(n+1):dim(x)[3],] <- a
        return(c)
      }
    }
    
    else{
      warning('Not implemented')
    }
  }
  
}

