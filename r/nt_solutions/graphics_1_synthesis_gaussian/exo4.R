




for (plt in 1:4){
  w <- array(rnorm(n*n), c(n,n)) / n
  w <- w - mean(w) + 1/n**2
  
  
  w_transf <- fft(w)
  w_transf_extended <- array(0, c(n,n,1,dim(f0)[4]))
  for (k in 1:dim(f0)[4]){
    w_transf_extended[,,,k] <- w_transf
  }
  
  FFT <- apply(f0, c(3,4), fft)
  CONV <- array(0, dim(f0))
  for (i in 1:n){ for (j in 1:n){ CONV[i,j,,] <- FFT[(i-1)*n + j,,] } }
  
  CONV <- CONV*w_transf_extended
  
  
  IFFT <- apply(CONV, c(3,4), fft, inverse=TRUE) / (n*n)
  f <- array(0, dim(f0))
  for (i in 1:n){ for (j in 1:n){ f[i,j,,] <- Re(IFFT[(i-1)*n + j,,]) } }
  
  u <- f
  u[1,1,,] <- c(0,0,0)
  u[2,1,,] <- c(1,1,1)
  imageplot(as.cimg(clamp(u)), 'Synthesized', c(2,2,plt))
}