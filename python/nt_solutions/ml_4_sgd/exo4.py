tau = .002/n
ElistG = zeros( ( int(niter/err_rate), nsamples) )
for irun in arange(0,nsamples):
    w = zeros( (p,1) )
    G = zeros( (p,n) ) # keep track of gradients
    g = zeros( (p,1) )
    for it in arange(0,niter):
        if mod(it,err_rate)==1:
            ElistG[ int(it/err_rate),irun ] = E(w,X,y)
        i = int( floor(random.rand()*n) )   # draw uniformly
        g1 = nablaEi(w,i)
        # update grad
        g = g - MakeCol(G[:,i]) + g1
        G[:,i] = g1.flatten()
        #
        w = w - tau * g

vmin = np.min( (np.min(Elist), ElistS.flatten().min(), ElistA.flatten().min(), ElistG.flatten().min() ) )

u = log10(ElistS-vmin+1e-20)
v = log10(ElistA -vmin+1e-20)
w = log10(ElistG -vmin+1e-20)

clf
plot(1,Inf, 'b')
plot(1,Inf, 'r')
plot(1,Inf, 'g')
plot( arange(0,niter,err_rate), u, 'b' )
plot( arange(0,niter,err_rate), v, 'r' )
plot( arange(0,niter,err_rate), w, 'g' )
axis((1,niter, np.min(w), np.max(w) ))
title('$log(E(w_l) - min E)$')
legend( ('SGD', 'SGA', 'SAG') )
