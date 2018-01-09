


cR <- as.vector(abs(fW))[order(-as.vector(abs(fW)))]
err_wav <- norm(f)**2 - cumsum(cR**2)
for (i in 1:length(err_wav)){ err_wav[i] <- max(err_wav[i], 1e-10)}

plot(log10(err_fft/norm(f)**2), type="l", col="red", lwd=2, xlim=c(1,n**2/50), ylim=c(-2.35,0))
lines(log10(err_wav/norm(f)**2), type="l", col="blue", lwd=2)
title(main="log_10(epsilon^2[M]/ ||f||^2)")
legend(x="topright",
       legend=c("Fourier", "Wavelets"),
       col=c("red", "blue"),
       inset = 0.05,
       lwd=2,
       cex=1)
