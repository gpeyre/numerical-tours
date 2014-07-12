niter = 600;
z = randn(n,n,2);
r = b-A(z);
p = r;
C = []; F = [];
for i=1:niter
    %
    x = z(:,:,1);
    C(i) = norm(y-Phi(x));
    F(i) = norm(delta(x));
    %
    r1 = r;
    Ap = A(p);
    alpha = dotp(r,r) / dotp(p,Ap);
    z = z + alpha*p;
    r = r-alpha*Ap;
    beta = dotp(r,r)/dotp(r1,r1);
    p = r + beta*p;
end
clf;
subplot(2,1,1);
plot(log10(C));
title('log_{10}|C(x_k)|');
axis tight;
subplot(2,1,2);
plot(log(abs(F(1:end/2)-F(end))));
title('log_{10}|F(x_k) - F(x^{*})|');
axis tight;
