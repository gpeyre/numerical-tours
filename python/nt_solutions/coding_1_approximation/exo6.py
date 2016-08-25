plt.figure(figsize = (8,8))
Mlist = [round(n**2/100), round(n**2/20)]

for i in range(len(Mlist)):
    M = Mlist[i]
    #perform hard thresholding 
    a = np.sort(np.ravel(abs(fC)))[::-1] #sort a 1D copy of F in descending order
    T = a[M]
    fCT = np.multiply(fC,abs(fC) > T)
    fM = idct2(fCT)
    imageplot(clamp(fM), "M/ N = %.2f, SNR = %.1f dB" %(M/n**2,snr(f,fM)),[1, 2, i+1])