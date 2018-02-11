tlist = np.linspace(.5, 2, 30)*sigma

snr_stein = [snr(f0, ThreshWav(f,t)) for t in Tlist]

plt.plot([T/sigma for T in Tlist], snr_stein, linewidth = 2)
plt.xlabel('$T/\sigma$')
plt.ylabel("SNR")

plt.show()

Tmax = np.max(snr_stein)
fBlock = ThreshWav(f, T)