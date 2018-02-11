



for (u in 1:length(Mlist)){
  M <- Mlist[u]
  a <- as.vector(abs(fL))[order(-as.vector(abs(fL)))] #sort a 1D copy of F in descending order
  T <- a[M]
  fLT <- fL * (abs(fL) > T)
  #fLT = perform_thresholding(fL, M, 'largest')
  fM <- fLT
  
  for (i in 1:(n%/%w)){
    for (j in 1:(n%/%w)){
      fM[((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w)] <- idct2(fLT[((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w)])
    }
  }
  
  imageplot(clamp(fM), paste("M/N =", format(M/n**2, digits=3),", SNR =", format(snr(f,fM), digits=3),"dB"), c(1, 2, u))
}