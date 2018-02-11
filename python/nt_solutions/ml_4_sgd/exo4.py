tau = .002/n
ElistG = np.zeros( ( int(niter/err_rate), nsamples) )
for irun in np.arange(0,nsamples):
    w = np.zeros( (p,1) )
    G = np.zeros( (p,n) ) # keep track of gradients
    g = np.zeros( (p,1) )
    for it in np.arange(0,niter):
        if np.mod(it,err_rate)==1:
            ElistG[ int(it/err_rate),irun ] = E(w,X,y)
        i = int( np.floor(np.random.rand()*n) )   # draw uniformly
        g1 = nablaEi(w,i)
        # update grad
        g = g - MakeCol(G[:,i]) + g1
        G[:,i] = g1.flatten()
        #
        w = w - tau * g

vmin = np.min( (np.min(Elist), ElistS.flatten().min(), ElistA.flatten().min(), ElistG.flatten().min() ) )

u = np.log10(ElistS-vmin+1e-20)
v = np.log10(ElistA -vmin+1e-20)
w = np.log10(ElistG -vmin+1e-20)

plt.clf
plt.plot(1,np.Inf, 'b')
plt.plot(1,np.Inf, 'r')
plt.plot(1,np.Inf, 'g')
plt.plot( np.arange(0,niter,err_rate), u, 'b' )
plt.plot( np.arange(0,niter,err_rate), v, 'r' )
plt.plot( np.arange(0,niter,err_rate), w, 'g' )
plt.axis((1,niter, np.min(w), np.max(w) ))
plt.title('$log(E(w_l) - min E)$')
plt.legend( ('SGD', 'SGA', 'SAG') )
