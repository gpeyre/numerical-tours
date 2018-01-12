niter = 100000
nsamples = 5
err_rate = 50

ElistS = zeros( ( int(niter/err_rate), nsamples) )
for irun in arange(0,nsamples):
    w = zeros( (p,1) )
    for it in arange(0,niter):
        if mod(it,err_rate)==1:
            ElistS[ int(it/err_rate),irun ] = E(w,X,y)
        tau = tau0 / (1+it/l0)
        i = int( floor(random.rand()*n) )   # draw uniformly
        w = w - tau * nablaEi(w,i)


subplot(2,1,1)
plot( arange(0,niter,err_rate), ElistS, 'b')
plot(1+arange(0,Elist.size)*n, Elist, 'k--')
axis((1,niter, ElistS.flatten().min(), ElistS.flatten().max() ))
title('$E(w_l)$')

vmin = np.min( (np.min(Elist), ElistS.flatten().min() ) )
u = log10(ElistS-vmin+1e-20)
v = log10(Elist -vmin+1e-20)
subplot(2,1,2)
plot( arange(0,niter,err_rate), u, 'b' )
plot(1+arange(0,Elist.size)*n,  v, 'k--' )
axis((1,niter, np.min(u), np.max(u) ))
title('$log(E(w_l) - min E)$')
