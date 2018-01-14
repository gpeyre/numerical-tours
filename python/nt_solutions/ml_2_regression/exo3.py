niter = 400
Jlist = zeros((niter,1))
w = zeros((p,1))
for i in arange(0,niter):
    Jlist[i] = J(w,Lambda)
    w = ISTA(w,Lambda,tau)
ndisp = int(niter/4)

clf
subplot(2,1,1)
plot(Jlist[0:ndisp])
axis('tight')
title('J(w_k)')
subplot(2,1,2)
e = log10( Jlist[0:ndisp]-Jlist.min() +1e-20)
plot(e-e[0])
axis('tight')
title('$log(J(w_k)-min J)$')
