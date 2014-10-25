gamma = gamma0;
displist = around(linspace(0,niter,10))
k = 0;
clf; 
imageplot(transpose(W));
for i in arange(0,niter+1):
    N = normal(gamma);
    g = EvalW(gamma) * normalC(gamma) - dotp(EvalG(gamma), N) * N;
    gamma = resample( gamma + dt*g );    
    if i==displist[k]:   
        lw = 1;    
        if i==0 or i==niter:
            lw = 4;
        cplot(gamma, 'r', lw);
        k = k+1;
        axis('equal'); axis('off');