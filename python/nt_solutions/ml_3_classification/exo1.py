niter = 5000
Elist = np.zeros([niter,1])
w = np.zeros( (p+1,1) ) # initialization
for i in np.arange(0,niter):
    w = w - tau * nablaE(w,AddBias(X),y)
    Elist[i] = E(w,AddBias(X),y)
ndisp = int(niter/2)

## display ##
plt.clf
plt.subplot(2,1,1)
plt.plot(np.arange(0,ndisp), Elist[0:ndisp])
plt.axis('tight')
plt.title('$E(w_l)$')
plt.subplot(2,1,2);
plt.plot(np.arange(10,ndisp), np.log10(Elist[10:ndisp]-np.min(Elist))  )
plt.axis('tight')
plt.title('$log(E(w_l) - min E)$');
