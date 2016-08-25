plt.figure(figsize = (8,8))
Mlist = [round(n**2/100), round(n**2/20)]

for i in range(len(Mlist)):
    M = Mlist[i]
    a = np.sort(np.ravel(abs(fW)))[::-1] #sort a 1D copy of F in descending order
    T = a[M]
    fWT = np.multiply(fW,(abs(fW) > T))
    fM = np.real(perform_wavortho_transf(fWT, Jmin, -1, h))
    imageplot(clamp(fM), "M/N = %.2f, SNR = %.1f dB" %(M/n**2,snr(f,fM)),[1, 2, i+1])