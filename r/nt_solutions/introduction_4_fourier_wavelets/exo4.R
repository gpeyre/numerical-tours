fw = perform_wavelet_transf(f, Jmin + 1, +1)
a = sort(Mod(F), decreasing=TRUE) #sort a 1D copy of F in descending order
T = a[M]
FT = F * (abs(F) > Re(T))
fM = perform_wavelet_transf(fw1, Jmin + 1, -1)
# Normalize to 0-1
fM = (fM-min(fM))/(max(fM)-min(fM))

imageplot(fM, paste("Non-linear, Wavelets, SNR = ",  round(snr(f, fM), 1), "dB"))


options(repr.plot.width=5, repr.plot.height=4)
plot(f[, n/2], type="l", col=4, ylab="", xlab="", main="f")
plot(fM[, n/2], type="l", col=4, ylab="", xlab="", main="fM")