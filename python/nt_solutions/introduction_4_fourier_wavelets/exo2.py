F = pyl.fft2(f)
a = np.sort(np.ravel(abs(F)))[::-1] #sort a 1D copy of F in descending order
T = a[M]
FT = np.multiply(F,(abs(F) > T))
fM = np.real(pyl.ifft2(FT))
plt.figure(figsize = (5,5))
imageplot(clamp(fM), "Non-linear, Fourier, SNR = %.1f dB" %snr(f, fM))