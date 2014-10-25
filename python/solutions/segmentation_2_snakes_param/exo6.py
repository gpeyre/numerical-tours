G = grad(W);
G = G[:,:,0] + 1j*G[:,:,1]
EvalG = lambda gamma: bilinear_interpolate(G, imag(gamma), real(gamma))
EvalW = lambda gamma: bilinear_interpolate(W, imag(gamma), real(gamma))
#
gamma = gamma0;
displist = around(linspace(0,niter,10))
k = 0;
clf
imageplot(transpose(f))
for i in arange(0,niter+1):
    N = normal(gamma)
    g = EvalW(gamma) * normalC(gamma) - dotp(EvalG(gamma), N) * N
    gamma = gamma + dt*g
    gamma = resample( gamma )
    # impose start/end point
    gamma[0] = x0
    gamma[-1] = x1
    if i==displist[k]:   
        lw = 1;    
        if i==0 or i==niter:
            lw = 4;
        cplot(gamma, 'r', lw);
        k = k+1;
        axis('equal'); axis('off');
        plot(real(gamma[0]), imag(gamma[0]), 'b.', markersize=20)
        plot(real(gamma[-1]), imag(gamma[-1]), 'b.', markersize=20);