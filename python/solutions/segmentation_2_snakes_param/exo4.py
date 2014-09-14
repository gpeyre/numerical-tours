G = grad(W);
G = G[:,:,0] + 1j*G[:,:,1]
EvalG = lambda gamma: bilinear_interpolate(G, imag(gamma), real(gamma))
EvalW = lambda gamma: bilinear_interpolate(W, imag(gamma), real(gamma))
#
gamma = gamma0
displist = around(linspace(0,niter,10))
k = 0;
clf
imageplot(transpose(f))
for i in arange(0,niter+1):
    n = normal(gamma)
    g = EvalW(gamma) * normalC(gamma) - dotp(EvalG(gamma), n) * n
    gamma = resample( gamma + dt*g )
    if i==displist[k]:   
        lw = 1;    
        if i==0 or i==niter:
            lw = 4;
        cplot(gamma, 'r', lw);
        k = k+1;
        axis('equal'); axis('off');