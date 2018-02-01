




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

