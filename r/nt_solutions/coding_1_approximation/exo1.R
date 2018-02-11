



Mlist <- c(round(n**2/100), round(n**2/20))

for (i in 1:length(Mlist)){
  M <- Mlist[i]
  #perform hard thresholding
  F <- fft(f)
  a <- as.vector(abs(F))[order(-as.vector(abs(F)))] #sort a 1D copy of F in descending order
  T <- a[M]
  fFT <- F * (abs(F) > T)
  fM <- Re(fft(fFT, inverse=TRUE)/length(fFT))
  imageplot(clamp(fM), paste("M/N =", format(M/n**2, digits=3),", SNR =", format(snr(f,fM), digits=3),"dB"), c(1, 2, i))
}