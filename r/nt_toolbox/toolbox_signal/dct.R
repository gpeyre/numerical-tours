



dct <- function (x, inverted=FALSE) 
{
  ####
  # Compute the dct/idct of vector/matrix x
  ####
  
  
  
  n <- length(x)
  res <- x
  
  if (!inverted){
    for (k in 0:(n-1)) {
      res[k+1] <- 2*sum(x*cos(pi/n*((0:(n-1))+0.5)*k))
    }
    res <- sqrt(1/(2*n))*res
    res[1] <- sqrt(1/2)*res[1]
  }
  
  if (inverted) {
    for (k in 0:(n-1)) {
      res[k+1] <- x[1]/sqrt(n) + sqrt(2/n)*sum(x[2:n]*cos(pi*(k+0.5)*(1:(n-1))/n))
    }
  }
  
  
  return(res)
  
}




dct_2d <- function(x, inverted=FALSE){
  
  return(t(apply(x, 1, dct, inverted=inverted)))
  
}
