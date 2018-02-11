

ntests <- 20
slist <- seq(0.01, 1.5, length.out=20)
err <- c()

for (i in 1:ntests){
  h <- exp(-(X**2 + Y**2 + Z**2)/(2*(slist[i]**2)))
  h <- h/sum(h)
  Mh <- Re(fft( fft(Mnoisy) * fft(fftshift_3d(h)), inverse=T) / length(h))
  #Mh <- rescale(Mh)
  err <- c(err, snr(M, Mh))
  if ( (i>2) & (err[i]>max(err[1:(i-1)])) ){
    Mblur <- Mh
  }
}


plot(slist, err, 'b', xlab="s", ylab="SNR")