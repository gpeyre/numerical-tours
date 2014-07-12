niter = 300;
x = zeros(n,1);
r = b-A*x;
p = r;
E = [];
for i=1:niter
    r1 = r;
    E(i) = norm(r(:));
    Ap = A*p;
    alpha = dotp(r,r) / dotp(p,Ap);
    x = x + alpha*p;
    r = r-alpha*Ap;
    beta = dotp(r,r)/dotp(r1,r1);
    p = r + beta*p;
end
clf;
plot(log10(E));
title('log_{10}|r_k|');
axis tight;
