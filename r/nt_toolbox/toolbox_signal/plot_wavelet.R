




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

