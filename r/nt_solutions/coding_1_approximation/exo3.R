

err_fft <- norm(f)**2 - cumsum(cR**2)
for (i in 1:length(err_fft)){ err_fft[i] <- max(err_fft[i], 1e-10)}

plot(log10(err_fft/norm(f)**2), type="l", col="blue", lwd=2, xlim=c(1,n**2/50), ylim=c(-2.35,0))
title(main="log_10(epsilon^2[M]/ ||f||^2)")


