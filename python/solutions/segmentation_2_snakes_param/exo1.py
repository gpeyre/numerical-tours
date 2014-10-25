gamma = gamma1
displist = around(linspace(0,niter,10))
k = 0;
for i in arange(0,niter+1):
    gamma = resample( gamma + dt * normalC(gamma) );
    if i==displist[k]:
        lw = 1;
        if i==0 or i==niter:
            lw = 4;
        cplot(gamma, 'r', lw);
        k = k+1;
        axis('tight');  axis('off');