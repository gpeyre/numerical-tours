F = fft(f)
a = sort(Mod(F), decreasing=TRUE) #sort a 1D copy of F in descending order
T = a[M]
FT = F * (abs(F) > Re(T))
fM = Re(fft(FT, inverse=TRUE))
# Normalize to 0-1
fM = (fM-min(fM))/(max(fM)-min(fM))

imageplot(fM, paste("Non - Linear, Fourier, SNR = ",  round(snr(f, fM), 1), "dB"))