


perform_thresholding <- function(f,M,type){
  ####
  # Only 3 types of thresholding currently implemented
  ####
  
  if (type == "largest"){
    a <- abs(as.vector(MW))
    a <- a[order(-a)]
    T <- a[M+1]
    y <- f * (abs(f) > T)
  }
  else if (type == "soft"){
    s <- abs(f) - M
    s <- (s + abs(s))/2
    y <- sign(f) * s
  }
  else if (type == "hard"){
    y <- f * (abs(f) > M)
  }
  return(y)
}