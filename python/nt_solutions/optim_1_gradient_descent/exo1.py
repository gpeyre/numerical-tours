x = x0
niter = 20
E = zeros((niter,1))
X = zeros((2,niter))
for i in arange(0,niter):
    X[:,i] = x.flatten()
    E[i] = f(x)
    x = x - tau*Gradf(x)

clf
h = plot(log10(E))
axis('tight')
title('log_{10}(x^{(k)})')
