

Tlist <- seq(0.8, 4.5, length=25)*sigma
err_soft <- matrix(0, c(length(Tlist),1))
err_hard <- matrix(0, c(length(Tlist),1))
for (i in 1:length(Tlist)){
  aT <- matrix(sapply(a, thresh_hard, t=Tlist[i]), c(n,n))
  fWav <- perform_wavortho_transf(aT,Jmin,-1,h)
  err_hard[i] <- snr(f0,fWav)
  aT <- matrix(sapply(a, thresh_soft, t=Tlist[i]), c(n,n))
  aT[seq(1,n,2^Jmin),seq(1,n,2^Jmin)] <- a[seq(1,n,2^Jmin),seq(1,n,2^Jmin)]
  fWav <- perform_wavortho_transf(aT,Jmin,-1,h)
  err_soft[i] <- snr(f0,fWav)
}
plot(Tlist/sigma,err_soft, type="l", col="green") 
lines(Tlist/sigma,err_hard, col="blue")
legend("topright",legend=c("Soft", "Hard"), col=c("green", "blue"), lty=c(1,1), cex=1)