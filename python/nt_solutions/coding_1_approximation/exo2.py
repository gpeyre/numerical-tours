cR = np.sort(np.ravel(abs(fF)))[::-1]
h = plt.plot(np.log10(cR), linewidth=2)
plt.xlim(0,n**2)

plt.show()