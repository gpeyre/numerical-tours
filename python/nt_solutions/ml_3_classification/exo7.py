niter = 1000
W = zeros((p,k))
Elist = zeros((niter,1))
tau = .01
for i in arange(0,niter):
    W = W - tau * nablaE(W)
    Elist[i] = E(W)
## display ##
clf
plot(arange(0,niter), Elist)
axis('tight')
