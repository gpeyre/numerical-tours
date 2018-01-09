



cR <- as.vector(abs(fL))[order(-as.vector(abs(fL)))]
err_ldct <- norm(f)**2 - cumsum(cR**2)
for (i in 1:length(err_ldct)){ err_ldct[i] <- max(err_ldct[i], 1e-10)}

plot(log10(err_fft/norm(f)**2), type="l", col="red", lwd=2, xlim=c(1,n**2/50), ylim=c(-2.35,0))
lines(log10(err_wav/norm(f)**2), type="l", col="blue", lwd=2)
lines(log10(err_dct/norm(f)**2), type="l", col="purple", lwd=2)
lines(log10(err_ldct/norm(f)**2), type="l", col="orange", lwd=2)

title(main="log_10(epsilon^2[M]/ ||f||^2)")
legend(x="topright",
       legend=c("Fourier", "Wavelets", "DCT", "Local DCT"),
       col=c("red", "blue", "purple", "orange"),
       inset = 0.05,
       lwd=2,
       cex=1)
