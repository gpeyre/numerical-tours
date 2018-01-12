tau0 = .05
ell0 = 100

#nsamples = 10
#err_rate = 50

ElistA = zeros( ( int(niter/err_rate), nsamples) )
for irun in arange(0,nsamples):
    w = zeros( (p,1) )
    w1 = w
    for it in arange(0,niter):
        if mod(it,err_rate)==1:
            ElistA[ int(it/err_rate),irun ] = E(w,X,y)
        tau = tau0 / (1+sqrt(it/ell0))
        i = int( floor(random.rand()*n) )   # draw uniformly
        w1 = w1 - tau * nablaEi(w1,i)
        w = 1.0/(it+1)*w1 + (1-1.0/(it+1))*w

# display
vmin = np.min( (np.min(Elist), ElistS.flatten().min(), ElistA.flatten().min() ) )
u = log10(ElistS-vmin+1e-20)
v = log10(ElistA -vmin+1e-20)
clf
plot(1,Inf, 'b')
plot(1,Inf, 'r')
plot( arange(0,niter,err_rate), u, 'b' )
plot( arange(0,niter,err_rate), v, 'r' )
axis((1,niter, np.min(v), np.max(v) ))
title('$log(E(w_l) - min E)$')
legend( ("SGD", "SGA") )
