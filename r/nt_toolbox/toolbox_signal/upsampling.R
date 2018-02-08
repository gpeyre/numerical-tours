



upsampling <- function(x, d){
  ####
  # up-sampling along dimension d by factor p=2
  ####
  p <- 2
  s <- dim(x)
  if (d==1){
    y <- matrix(0, p*s[1], s[2])
    y[ ((1:dim(y)[1])%%p == 1 ) , ] <- x
  }
  else if (d==2){
    y <- matrix(0, s[1], p*s[2])
    y[ , ((1:dim(y)[2])%%p == 1 ) ] <- x
  }
  else{
    warning('Not implemented')
  }
  return(y)
}

