




gaussian_blur <- function(f, sigma){
  ####
  # gaussian_blur - gaussian blurs an image
  #
  # M = perform_blurring(M, sigma, options);
  #
  # M is the original data
  # sigma is the std of the Gaussian blur (in pixels)
  # 
  # Copyright (c) 2007 Gabriel Peyre
  ####
  if (sigma<=0){ return() }
  n <- max(dim(f))
  t <- c(0:(n/2),(-n/2):(-2))
  Y <- meshgrid_2d(t,t)$X ; X <- meshgrid_2d(t,t)$Y
  h <- exp(-(X**2 + Y**2)/(2.0*sigma**2))
  h <- h/sum(h)
  fft_prod <- fft(f)*fft(h)
  return( Re( fft(fft_prod, inverse=T)/length(fft_prod) ) )
  
}

