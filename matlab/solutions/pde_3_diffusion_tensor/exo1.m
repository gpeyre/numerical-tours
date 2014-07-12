t = 15;
niter=  round(t/tau);
kdisp = round(linspace(0,niter,5)); kdisp(1) = [];
k =1;
f1 = f;
for i=1:niter
    f1 = f1 + tau * div( Mult(S, nabla(f1) ) );
    if i==kdisp(k)
        subplot(2,2,k);
        imageplot(clamp(f1));
        k = k+1;
    end
end
