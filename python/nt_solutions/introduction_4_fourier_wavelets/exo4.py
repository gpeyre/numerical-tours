fw = perform_wavelet_transf(f, Jmin, +1)
a = np.sort(np.ravel(abs(fw)))[::-1] #sort a 1D copy of fw in descending order
T = a[M]
fw1 = np.multiply(fw,(abs(fw) > T))
fM = perform_wavelet_transf(fw1, Jmin, -1)
plt.figure(figsize=(5,5))
imageplot(clamp(fM), "Non-linear, Wavelets, SNR = %.1f dB" %snr(f, fM))