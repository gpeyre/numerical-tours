Pi = @(f)f.*(1-Gamma) + y.*Gamma;
tau = .8/4;
niter = 300;
E = [];
k = 1; ndisp = [1 5 10 niter];
norm1 = @(f)norm(f(:));
f = y; f(Gamma==0) = .5;
for i=1:niter
    E(i) = norm1(grad(f));
    f = Pi( f + tau*Delta(f) );
end
clf;
imageplot(f);
