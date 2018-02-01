



read_bin <- function(name, ndims = 2){
  ####
  # reading from a binary file in 3 dimensions
  ####
  file_dim <- readBin(name, what="raw", n=2*ndims)
  
  if (ndims == 2){
    n <- readBin(file_dim[1], what="integer", size=1, n=1, signed=FALSE)
    p <- readBin(file_dim[3], what="integer", size=1, n=1, signed=FALSE)
    file_raw <- readBin(name, what="raw", n=n*p + 2*ndims)
    file_trans <- readBin(file_raw[5:length(file_raw)], what="integer", size=1, n=n*p, signed=FALSE)
    dim(file_trans) <- c(n,p)
  }
  else{
    n <- readBin(file_dim[1], what="integer", size=1, n=1, signed=FALSE)
    p <- readBin(file_dim[3], what="integer", size=1, n=1, signed=FALSE)
    q <- readBin(file_dim[5], what="integer", size=1, n=1, signed=FALSE)
    file_raw <- readBin(name, what="raw", n=n*p*q + 2*ndims)
    file_trans <- readBin(file_raw[7:length(file_raw)], what="integer", size=1, n=n*p*q, signed=FALSE)
    dim(file_trans) <- c(n,p,q)
  }
  return(file_trans)
}