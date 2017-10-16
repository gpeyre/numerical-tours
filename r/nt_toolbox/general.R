



# Functions


norm <- function(v){
  ####
  # Euclidian norm of a vector
  ####
  return(sqrt(sum(v**2)))
}



roll <- function(x, n){
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
    else{
      warning('Not implemented')
    }
  }
  
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



clamp <- function(x, a=c(), b=c()){
  ####
  # clamp - clamp a value
  # 
  # y = clamp(x,a,b);
  # 
  # Default is [a,b]=[0,1].
  # 
  # Copyright (c) 2004 Gabriel Peyre
  ####
  if (length(a)==0){
    a <- 0.0
  }
  if (length(b)==0){
    b <- 1.0
  }
  return( pmin( pmax(x,a), b) )
}




ravel <- function(M){
  ####
  # return the 1D-array correponding to matrix M (row-wise)
  ####
  return( as.vector(t(M)) )
}





