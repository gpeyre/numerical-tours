from pylab import linalg

err_fft = [max(e,1e-10) for e in linalg.norm(f)**2 - np.cumsum(cR**2)]
plt.plot(np.log10(err_fft/linalg.norm(f)**2),linewidth=2)
plt.title("$\log_{10}(\epsilon^2[M]/ ||f||^2)$")
plt.xlim(1,n**2/50)
plt.ylim(-2.35,0)

plt.show()