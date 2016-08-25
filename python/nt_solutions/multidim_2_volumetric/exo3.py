ntests = 20
slist = np.linspace(.01, 1.5, ntests)
err = []
for i in range(ntests):
    h = np.exp(-(X**2 + Y**2 + Z**2)/(2*slist[i]**2))
    h = h/np.sum(h)
    Mh = np.real(pyl.ifft2(pyl.fft2(Mnoisy,axes=(0,1,2)) * pyl.fft2(pyl.fftshift(h,axes= (0,1,2)),axes=(0,1,2)),axes= (0,1,2)))
    err = err + [snr(M, Mh)]
    if i > 1 and err[i] > np.max(err[:i]):
        Mblur = Mh

plt.figure(figsize=(7,5))
plt.plot(slist, err, '.-')
plt.xlabel("s")
plt.ylabel("SNR")
plt.show()