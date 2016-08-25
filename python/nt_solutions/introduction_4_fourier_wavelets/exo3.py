q = int(np.sqrt(M))
fw = perform_wavelet_transf(f, Jmin, +1)
fw1 = np.zeros([n0, n0])
fw1[:q, :q] = fw[:q, :q]
fM = perform_wavelet_transf(fw1, Jmin, -1)
plt.figure(figsize=(5,5))
imageplot(clamp(fM),  "Linear, Wavelets, SNR = %.1f dB" %snr(f, fM))