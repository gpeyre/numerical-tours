niter = 5000
Elist = zeros([niter,1])
w = zeros( (p+1,1) ) # initialization
for i in range(0,niter):
    w = w - tau * nablaE(w,AddBias(X),y)
    Elist[i] = E(w,AddBias(X),y)
ndisp = int(niter/2)

## display ##
clf
subplot(2,1,1)
plot(arange(0,ndisp), Elist[0:ndisp])
axis('tight')
title('$E(w_l)$')
subplot(2,1,2);
plot(arange(1,ndisp), log10(Elist[1:ndisp]-np.min(Elist))  )
axis('tight')
title('$log(E(w_l) - min E)$');
