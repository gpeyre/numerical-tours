niter = 400
flist = np.zeros((niter,1))
x = np.zeros((p,1))
for i in np.arange(0,niter):
    flist[i] = f(x,Lambda)
    x = ISTA(x,Lambda,tau)
ndisp = int(niter/4)

plt.clf
plt.subplot(2,1,1)
plt.plot(flist[0:ndisp])
plt.axis('tight')
plt.title('f(x_k)')
plt.subplot(2,1,2)
e = np.log10( flist[0:ndisp]-flist.min() +1e-20)
plt.plot(e-e[0])
plt.axis('tight')
plt.title('$log(f(x_k)-min f)$')
