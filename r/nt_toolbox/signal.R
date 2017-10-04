


# Libraries & Sources


dir_0 <- getwd()
setwd(dirname(parent.frame(2)$ofile))
source("general.R")

library(imager)






# Functions


cconv <- function(x, h, d){
  ####
  # Circular convolution along dimension d.
  # h should be small and with odd size
  ####
  if (d==2){
    # apply to transposed matrix
    return(t(cconv(t(x),h,1)))
  }
  y <- array(0, dim(x))
  p <- length(h)
  pc <- as.integer(round( as.numeric((p - 1) / 2 )))
  for (i in (1:p)){
    y <- y + h[i] * circshift1d(x, i - pc - 1)
  }
  return(as.matrix(y))
}



subsampling <- function(x, d){
  ####
  # subsampling along dimension d by factor p=2
  ####
  p <- 2
  if (d==1){
    y <- x[ ((1:dim(x)[1])%%p == p-1 ) , ]
  }
  else if (d==2){
    y <- x[ , ((1:dim(x)[2])%%p == p-1 ) ]
  }
  else{
    warning('Not implemented')
  }
  return(as.matrix(y))
}



upsampling <- function(x, d){
  ####
  # up-sampling along dimension d by factor p=2
  ####
  p <- 2
  s <- dim(x)
  if (d==1){
    y <- matrix(0, p*s[1], s[2])
    y[ ((1:dim(x)[1])%%p == p-1 ) , ] <- x
  }
  else if (d==2){
    y <- matrix(0, s[1], p*s[2])
    y[ , ((1:dim(x)[2])%%p == p-1 ) ] <- x
  }
  else{
    warning('Not implemented')
  }
  return(as.matrix(y))
}



load_image <- function(name, n=-1, flatten=1, resc=1, grayscale=1){
  ####
  # Load an image from a file, rescale its dynamic to [0,1],
  # turn it into a grayscale image and resize it to size n x n.
  ####
  f <- load.image(name)
  # turn into normalized grayscale image
  if (grayscale == 1){
    if ( (flatten==1) & (length(dim(f))>2) ){
      f <- apply(f, c(1,2), sum)
      f <- as.cimg(f)
    }
  }
  if (resc==1){
    f <- rescale(f)
  }
  # change the size of the image
  if (n>0){
    if (dim(f)[3]==1){
      f <- resize(f, size_x = n, size_y = n, interpolation_type = 5)
    }
    else if (dim(f)[3]>1){
      f <- resize(f, size_x = n, size_y = n, size_z = dim(f)[3], interpolation_type = 5)
    }
    f <- t( f[1:n,1:n] )
    return(as.cimg(f))
  }
}



imageplot <- function(f, str='', sbpt=c()){
  ####
  # Use nearest neighbor interpolation for the display.
  ####
  #f <- t( f[1:n,f:n] )
  f <- as.cimg( t(as.matrix(f)) )
  if (length(sbpt) >0){
    if (sbpt[3]==1){
      par(mfrow=sbpt[1:2]) 
    }
  }
  plot(f, interpolate = TRUE, colorscale = gray, axes = FALSE, main = str)
}





# Replacing the working directory

setwd(dir_0)