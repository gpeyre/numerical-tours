from pylab import linalg

cR = np.sort(np.ravel(abs(fW)))[::-1]
err_wav = [max(e,1e-10) for e in linalg.norm(f)**2 - np.cumsum(cR**2)]

plt.plot(np.log10(err_fft/linalg.norm(f)**2),linewidth=2, color = "red", label = "Fourier")
plt.plot(np.log10(err_wav/linalg.norm(f)**2),linewidth=2, color = "blue", label = "Wavelets")

plt.title("$\log_{10}(\epsilon^2[M]/ ||f||^2)$")
plt.xlim(1,n**2/50)
plt.ylim(-2.35,0)
plt.legend()

plt.show()