u = ones(N,1);
niter = 500;
gamma = .001;
pi = exp( -C/gamma );
E1 = []; E2 = [];
for i=1:niter
    pi = ProjC1(pi,p);
    E2(end+1) = norm( pi'*u-q )/norm(q);
    pi = ProjC2(pi,q);
    E1(end+1) = norm( pi*u-p )/norm(p);
end
% Display error decay.
clf;
plot( log10([E1; E2]') );
axis tight;
