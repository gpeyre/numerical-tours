



Mlist <- c(round(n**2/100), round(n**2/20))

for (i in 1:length(Mlist)){
  M <- Mlist[i]
  a <- as.vector(abs(fW))[order(-as.vector(abs(fW)))] #sort a 1D copy of F in descending order
  T <- a[M]
  fWT <- fW * (abs(fW) > T)
  fM <- Re(perform_wavortho_transf(fWT, Jmin, -1, h))
  imageplot(clamp(fM), paste("M/N =", format(M/n**2, digits=3),", SNR =", format(snr(f,fM), digits=3),"dB"), c(1, 2, i))
}