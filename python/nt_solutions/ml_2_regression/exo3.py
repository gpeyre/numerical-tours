niter = 400
Jlist = np.zeros((niter,1))
w = np.zeros((p,1))
for i in np.arange(0,niter):
    Jlist[i] = J(w,Lambda)
    w = ISTA(w,Lambda,tau)
ndisp = int(niter/4)

plt.clf
plt.subplot(2,1,1)
plt.plot(Jlist[0:ndisp])
plt.axis('tight')
plt.title('J(w_k)')
plt.subplot(2,1,2)
e = np.log10( Jlist[0:ndisp]-Jlist.min() +1e-20)
plt.plot(e-e[0])
plt.axis('tight')
plt.title('$log(J(w_k)-min J)$')
