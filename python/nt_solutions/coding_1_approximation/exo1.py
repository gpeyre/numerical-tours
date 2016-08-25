plt.figure(figsize = (8,8))
Mlist = [round(n**2/100), round(n**2/20)]

for i in range(len(Mlist)):
    M = Mlist[i]
    #perform hard thresholding
    F = pyl.fft2(f)
    a = np.sort(np.ravel(abs(F)))[::-1] #sort a 1D copy of F in descending order
    T = a[M]
    fFT = np.multiply(F,(abs(F) > T))
    fM = np.real(pyl.ifft2(fFT))
    imageplot(clamp(fM), "M/N = %.2f, SNR = %.1f dB" %(M/n**2,snr(f,fM)),[1, 2, i+1])