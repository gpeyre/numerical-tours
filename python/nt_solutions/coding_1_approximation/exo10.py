plt.figure(figsize = (8,8))

for u in range(len(Mlist)):
    M = Mlist[u]
    a = np.sort(np.ravel(abs(fL)))[::-1] #sort a 1D copy of F in descending order
    T = a[M]
    fLT = np.multiply(fL,abs(fL) > T)
    #fLT = perform_thresholding(fL, M, 'largest')
    fM = fLT
    
    for i in range(1,n//w+1):
        for j in range(1,n//w+1):
            fM[(i-1)*w: i*w, (j-1)*w: j*w] = idct2(fLT[(i-1)*w: i*w, (j-1)*w: j*w])
    
    imageplot(clamp(fM), "M/N = %.2f, SNR = %.1f dB" %(M/n**2,snr(f,fM)),[1, 2, u+1])