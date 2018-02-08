


subsampling <- function(x, d){
  ####
  # subsampling along dimension d by factor p=2
  ####
  p <- 2
  if (d==1){
    y <- x[ ((1:dim(x)[1])%%p == 1 ) , , drop=F]
  }
  else if (d==2){
    y <- x[ , ((1:dim(x)[2])%%p == 1 ) , drop=F]
  }
  else{
    warning('Not implemented')
  }
  return(y)
}


