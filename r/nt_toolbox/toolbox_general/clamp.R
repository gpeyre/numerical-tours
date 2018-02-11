




clamp <- function(x, a=0, b=1){
  ####
  # clamp - clamp a value
  #
  #   y = clamp(x,a,b);
  #
  # Default is [a,b]=[0,1].
  ####
  x[ x>b ] <- b
  x[ x<a ] <- a
  return(x)
}

