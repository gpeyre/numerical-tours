q = as.integer(sqrt(M))
F = fftshift(fft(f[,]))
Sel = matrix(0, n, n)

Sel[c((n / 2 - q / 2): (n / 2 + q / 2)), c((n / 2 - q / 2): (n / 2 + q / 2))] = 1

F_zeros = F * Sel

fM = Re(fft(fftshift(F_zeros), inverse=TRUE))
# Normalize to 0-1
fM = (fM-min(fM))/(max(fM)-min(fM))
imageplot(fM, paste("Linear, Fourier, SNR = ",  round(snr(f, fM), 1), "dB"))


options(repr.plot.width=5, repr.plot.height=4)
plot(f[, n/2], type="l", col=4, ylab="", xlab="", main="f")
plot(fM[, n/2], type="l", col=4, ylab="", xlab="", main="fM")
