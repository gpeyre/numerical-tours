niter = 2000
Flist = np.zeros([niter,1])
tau = .5
h = np.zeros([n,1])
for i in np.arange(0,niter):
    h = h - tau * nablaF(h, K, y)
    Flist[i] = F(h, K, y)
plt.clf
plt.plot(np.arange(0,niter), Flist)
plt.axis('tight')
