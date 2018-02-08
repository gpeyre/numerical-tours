



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

