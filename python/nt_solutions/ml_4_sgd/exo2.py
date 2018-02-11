niter = 100000
nsamples = 5
err_rate = 50

ElistS = np.zeros( ( int(niter/err_rate), nsamples) )
for irun in np.arange(0,nsamples):
    w = np.zeros( (p,1) )
    for it in np.arange(0,niter):
        if np.mod(it,err_rate)==1:
            ElistS[ int(it/err_rate),irun ] = E(w,X,y)
        tau = tau0 / (1+it/l0)
        i = int( np.floor(np.random.rand()*n) )   # draw uniformly
        w = w - tau * nablaEi(w,i)


plt.subplot(2,1,1)
plt.plot( np.arange(0,niter,err_rate), ElistS, 'b')
plt.plot(1+np.arange(0,Elist.size)*n, Elist, 'k--')
plt.axis((1,niter, ElistS.flatten().min(), ElistS.flatten().max() ))
plt.title('$E(w_l)$')

vmin = np.min( (np.min(Elist), ElistS.flatten().min() ) )
u = np.log10(ElistS-vmin+1e-20)
v = np.log10(Elist -vmin+1e-20)
plt.subplot(2,1,2)
plt.plot( np.arange(0,niter,err_rate), u, 'b' )
plt.plot(1+np.arange(0,Elist.size)*n,  v, 'k--' )
plt.axis((1,niter, np.min(u), np.max(u) ))
plt.title('$log(E(w_l) - min E)$')
