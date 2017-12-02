


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
  # if (n == -1){n = 512}
  if (n > 0){
    if (dim(f)[3]==1){
      f <- resize(f, size_x = n, size_y = n, interpolation_type = 5)
    }
    else if (dim(f)[4]>1){
      f <- resize(f, size_x = n, size_y = n, size_z = dim(f)[3], interpolation_type = 5)
      return(f)
    }
    f <- t(f[1:n,1:n])
    return(as.cimg(f))
  }
  return(as.cimg(f))

}



imageplot <- function(f, str='', sbpt=c()){
  ####
  # Use nearest neighbor interpolation for the display.
  ####
  #f <- t( f[1:n,f:n] )
  
  if ((length(dim(f)) > 3) && (dim(f)[4] > 1))
  {
        plot(f, main=str, axes=FALSE)
  }
  else
  {
  f <- as.cimg( t(as.matrix(f)) )
  if (length(sbpt) >0){
    if (sbpt[3]==1){
      par(mfrow=sbpt[1:2]) 
    }
  }
  plot(f, interpolate = FALSE, colorscale = gray, axes = FALSE, main = str)
  }
}





plot_wavelet <- function(fW, Jmin=0){
  ####
  # plot_wavelet - plot wavelets coefficients.
  
  # U = plot_wavelet(fW, Jmin):
    
  # Copyright (c) 2014 Gabriel Peyre
  ####
  
  
  rescaleWav <-function(A){
    v <- abs(max(A))
    B <- A
    if (v>0){
      B <- 0.5 + 0.5*(A/v)
    }
    return(B)
  }
  
  ##
  
  n <- dim(fW)[2]
  Jmax <- as.integer(log(n,2) - 1)
  U <- fW
  for ( j in Jmax:Jmin ){
    U[1:2**j, (2**j + 1):2**(j + 1)] <- rescaleWav(U[1:2**j, (2**j + 1):2**(j + 1)])
    U[(2 ** j + 1):2 ** (j + 1), 1:2**j] <- rescaleWav(U[(2 ** j + 1):2 ** (j + 1), 1:2**j])
    U[(2**j + 1):2**(j + 1), (2**j + 1):2**(j + 1)] <- (rescaleWav(U[(2**j + 1):2**(j + 1), (2**j + 1):2**(j + 1)]))
  }
  # coarse scale
  U[1:2**j, 1:2**j] <- rescale(U[1:2**j, 1:2**j])
  # plot underlying image
  imageplot(U)
  # display crosses
  for ( j in Jmax:Jmin ){
    points(c(1, 2**(j + 1)), c(2**j + 1, 2**j + 1), 'l', col="red")
    points(c(2**j, 2**j), c(0, 2**(j + 1)), 'l', col="red")
  }
  # display box
  points(c(1, n), c(1, 1), 'l', col="red")
  points(c(1, n), c(n, n), 'l', col="red")
  points(c(1, 1), c(1, n), 'l', col="red")
  points(c(n, n), c(1, n ), 'l', col="red")
  # return(U)
}





perform_wavortho_transf <- function(f, Jmin, dir, h){
  ####
  # perform_wavortho_transf - compute orthogonal wavelet transform
  
  # fw = perform_wavortho_transf(f,Jmin,dir,options);
  
  # You can give the filter in options.h.
  
  # Works in 2D only.
  
  # Copyright (c) 2014 Gabriel Peyre
  ####
  
  
  n <- dim(f)[2]
  Jmax <- as.integer(log(n, 2) - 1)
  # compute g filter
  u <- (- rep(1, length(h) - 1)) ** (1:(length(h)-1))
  # alternate +1/-1
  g <- c(0, rev(h[2:length(h)]) * u)
  
  if (dir == 1){
    ### FORWARD ###
    fW <- as.matrix(f)
    for (j in Jmax:Jmin){
      A <- fW[1:2**(j+1), 1:2**(j+1)]
      for (d in (1:2)){
        Coarse <- subsampling(cconv(A, h, d), d)
        Detail <- subsampling(cconv(A, g, d), d)
        if (d==1){
          A <- rbind(Coarse, Detail)
        }
        else{
          A <- cbind(Coarse, Detail)
        }
      }
      fW[1:2 ** (j + 1), 1:2 ** (j + 1)] <- A
    }
    return(fW)
  }
  else{
    ### BACKWARD ###
    fW <- f
    f1 <- fW
    for (j in Jmin:Jmax){
      A <- f1[1:2 ** (j + 1), 1:2 ** (j + 1)] ###
      for (d in (1:2)){
        if (d==1){
          Coarse <- A[1:2**j, , drop=F]
          Detail <- A[(2**j + 1): 2**(j + 1), , drop=F]
        }
        else{
          Coarse <- A[ , 1:2 ** j, drop=F]
          Detail <- A[ , (2 ** j + 1):2 ** (j + 1), drop=F]
        }
        Coarse <- cconv(upsampling(Coarse, d), reverse(h), d)
        Detail <- cconv(upsampling(Detail, d), reverse(g), d)
        A <- Coarse + Detail
      }
      f1[1:2 ** (j + 1), 1:2 ** (j + 1)] <- A
    }
    return(f1)
  }
}






snr <- function(x, y){
  ####
  # snr - signal to noise ratio
  # 
  #     v = snr(x,y);
  # 
  # v = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )
  # 
  #     x is the original clean signal (reference).
  #     y is the denoised signal.
  # 
  # Copyright (c) 2014 Gabriel Peyre
  ####
  x <- as.matrix(x)
  y <- as.matrix(y)
  return(20 * log( norm(x) / norm(x-y), 10) )
}


snr_2 <- function(x, y){
  ####
  # snr - signal to noise ratio
  # 
  #     v = snr(x,y);
  # 
  # v = 20*log10( norm(x(:)) / norm(x(:)-y(:)) )
  # 
  #     x is the original clean signal (reference).
  #     y is the denoised signal.
  # 
  # Copyright (c) 2014 Gabriel Peyre
  ####
  return(20 * log( norm(x) / norm(x-y), 10) )
}



div = function(x)
{ 
    # Divergence operator
    n = dim(x)[1]
    hdiff1 = x[c(2:n, 1), , 1]
    hdiff2 = x[, c(2:n, 1), 2]
    return(hdiff1 - x[,,1] + hdiff2 - x[,,2])
}

grad = function(x){
    n = dim(x)[1]
    hdiff = x[,] - x[c(n, 1:n-1),]
    vdiff = x[,] - x[,c(n, 1:n-1)]
    return (array(c(hdiff, vdiff), dim=c(n, n, 2)))
}


fftshift <- function(input_matrix, dim = -1) {

    rows <- dim(input_matrix)[1]    
    cols <- dim(input_matrix)[2]    

    swap_up_down <- function(input_matrix) {
        rows_half <- ceiling(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- ceiling(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    if (dim == -1) {
        input_matrix <- swap_up_down(input_matrix)
        return(swap_left_right(input_matrix))
    }
    else if (dim == 1) {
        return(swap_up_down(input_matrix))
    }
    else if (dim == 2) {
        return(swap_left_right(input_matrix))
    }
    else {
        stop("Invalid dimension parameter")
    }
}

ifftshift <- function(input_matrix, dim = -1) {

    rows <- dim(input_matrix)[1]    
    cols <- dim(input_matrix)[2]    

    swap_up_down <- function(input_matrix) {
        rows_half <- floor(rows/2)
        return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
    }

    swap_left_right <- function(input_matrix) {
        cols_half <- floor(cols/2)
        return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
    }

    if (dim == -1) {
        input_matrix <- swap_left_right(input_matrix)
        return(swap_up_down(input_matrix))
    }
    else if (dim == 1) {
        return(swap_up_down(input_matrix))
    }
    else if (dim == 2) {
        return(swap_left_right(input_matrix))
    }
    else {
        stop("Invalid dimension parameter")
    }
}







# Replacing the working directory

setwd(dir_0)