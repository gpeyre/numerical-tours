niter = 1000
W = np.zeros((p,k))
Elist = np.zeros((niter,1))
tau = .01
for i in np.arange(0,niter):
    W = W - tau * nablaE(W)
    Elist[i] = E(W)
## display ##
plt.clf
plt.plot(np.arange(0,niter), Elist)
plt.axis('tight')
