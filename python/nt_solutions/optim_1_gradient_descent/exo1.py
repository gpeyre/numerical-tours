x = x0
niter = 20
E = np.zeros((niter,1))
X = np.zeros((2,niter))
for i in np.arange(0,niter):
    X[:,i] = x
    E[i] = f(x)
    x = x - tau*Gradf(x)

plt.clf
h = plt.plot(np.log10(E))
plt.axis('tight')
plt.title('$log_{10}f(x^{(k)})$')
