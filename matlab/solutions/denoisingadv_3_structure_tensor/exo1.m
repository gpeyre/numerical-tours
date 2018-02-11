t = 20;
dt = .2;
niter=  round(t/dt);
kdisp = round(linspace(0,niter,5)); kdisp(1) = [];
k =1;
f1 = f; % rand(n);
for i=1:niter
    f1 = f1 + dt * div( TensorMult(T, grad(f1,options) ) );
    if i==kdisp(k)
        subplot(2,2,k);
        imageplot(f1);
        k = k+1;
    end
end
