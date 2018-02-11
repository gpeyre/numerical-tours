tau0 = .05
ell0 = 100

ElistA = np.zeros( ( int(niter/err_rate), nsamples) )
for irun in np.arange(0,nsamples):
    w = np.zeros( (p,1) )
    w1 = w
    for it in np.arange(0,niter):
        if np.mod(it,err_rate)==1:
            ElistA[ int(it/err_rate),irun ] = E(w,X,y)
        tau = tau0 / (1+np.sqrt(it/ell0))
        i = int( np.floor(np.random.rand()*n) )   # draw uniformly
        w1 = w1 - tau * nablaEi(w1,i)
        w = 1.0/(it+1)*w1 + (1-1.0/(it+1))*w

# display
vmin = np.min( (np.min(Elist), ElistS.flatten().min(), ElistA.flatten().min() ) )
u = np.log10(ElistS-vmin+1e-20)
v = np.log10(ElistA -vmin+1e-20)
plt.clf
plt.plot(1,np.Inf, 'b')
plt.plot(1,np.Inf, 'r')
plt.plot( np.arange(0,niter,err_rate), u, 'b' )
plt.plot( np.arange(0,niter,err_rate), v, 'r' )
plt.axis((1,niter, np.min(v), np.max(v) ))
plt.title('$log(E(w_l) - min E)$')
plt.legend( ("SGD", "SGA") )
