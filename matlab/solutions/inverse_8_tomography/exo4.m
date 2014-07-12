cb = @(x)reshape(Phi(PhiS(reshape(x,n,m))), n*m, 1);
niter = 10;
[y1,FLAG,RELRES,ITER,RESVEC] = cgs(cb,y(:),1e-10,niter);
y1 = reshape(y1,n,m);
fL2 = PhiS(y1);
clf;
imageplot(clamp(fL2));
