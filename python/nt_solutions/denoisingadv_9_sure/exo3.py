plt.figure(figsize = (5,5))
tmp = np.min(E)
i = np.argmin(E)
lambd = lambda_list[i]
hf = h(f,lambd)
imageplot(clamp(hf), "SNR = %.1f dB" %snr(f0, hf))