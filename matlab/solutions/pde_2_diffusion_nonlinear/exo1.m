lambda = 1e-2;
T = .5/lambda;
tau = .2;
niter = ceil(T/tau);
f = f0;
clf; k = 0;
for i=1:niter
    f = f + tau * div( g(A(f),lambda) .* grad(f) );
    if mod(i,floor(niter/4))==0
        k = k+1;
        imageplot(clamp(f), strcat(['T=' num2str(T*k/4,3)]), 2,2,k );
    end
end
