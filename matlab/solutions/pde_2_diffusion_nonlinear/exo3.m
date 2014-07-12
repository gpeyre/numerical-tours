T = 100;
tau = .2;
niter = ceil(T/tau);
f = f0;
clf; k = 0;
for i=1:niter
    f = f + tau * amplitude(grad(f)) .* curv(f);
    if mod(i,floor(niter/4))==0
        k = k+1;
        imageplot(clamp(f), strcat(['T=' num2str(T*k/4,3)]), 2,2,k );
    end
end
