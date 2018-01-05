niter = 2000
Flist = zeros([niter,1])
tau = .5
h = zeros([n,1])
for i in arange(0,niter):
    h = h - tau * nablaF(h, K, y)
    Flist[i] = F(h, K, y)
clf
plot(arange(0,niter), Flist)
axis('tight')
