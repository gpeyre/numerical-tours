
MW <- perform_haar_transf(Mnoisy, 1, +1)

Tlist <- seq(1, 4, length.out=20)*sigma
err_hard <- c()
err_soft <- c()

for ( i in (1:length(Tlist)) ){
  MWT <- perform_thresholding(MW, Tlist[i], 'hard')
  M1 <- perform_haar_transf(MWT, 1, -1)
  err_hard <- c(err_hard, snr(M, M1))
  MWT <- perform_thresholding(MW, Tlist[i], 'soft')
  M1 <- perform_haar_transf(MWT, 1, -1)
  err_soft <- c(err_soft, snr(M, M1))
  if ( (i>2) & (err_soft[i]>max(err_soft[1:(i-1)])) ){
    Mwav <- M1
  }
}

plot(Tlist/ sigma, err_hard, 'b', col="blue", cex=0.8,
     xlab="T/sigma", ylab="SNR", ylim=c(min(err_hard), max(err_soft)))
lines(Tlist/ sigma, err_soft, 'b', col="red", cex=0.8)
legend(x="topright",
       legend=c("hard", "soft"),
       col=c("blue", "red"),
       inset = 0.05,
       lwd=1,
       cex=0.8)