rho = .25;
niter = 12*4;
k = 0;
clf;
f1 = f;
for i=1:niter
    f1 = W(f1, rho*U );
    if mod(i,niter/4)==0
        k = k+1;
        imageplot(f1, strcat(['t=' num2str(i*rho)]), 2,2,k);
    end
end
