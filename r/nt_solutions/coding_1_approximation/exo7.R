



cR <- as.vector(abs(fC))[order(-as.vector(abs(fC)))]
err_dct <- norm(f)**2 - cumsum(cR**2)
for (i in 1:length(err_dct)){ err_dct[i] <- max(err_dct[i], 1e-10)}

plot(log10(err_fft/norm(f)**2), type="l", col="red", lwd=2, xlim=c(1,n**2/50), ylim=c(-2.35,0))
lines(log10(err_dct/norm(f)**2), type="l", col="green", lwd=2)
title(main="log_10(epsilon^2[M]/ ||f||^2)")
legend(x="topright",
       legend=c("Fourier", "DCT"),
       col=c("red", "green"),
       inset = 0.05,
       lwd=2,
       cex=1)
