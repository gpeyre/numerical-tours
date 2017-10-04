



# Functions


norm <- function(v){
  return(sqrt(sum(v**2)))
}



roll <- function(x, n){
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
  # Reverse a vecotr.
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


