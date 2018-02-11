% build vector field
x = linspace(-1,1,n);
[Y,X] = meshgrid(x,x);
d = sqrt(X.^2+Y.^2)+1e-5;
u = cat(3,-Y./d,X./d);
v = cat(3,X./d,Y./d);
% build tensor field
rho = 1/100;
T = perform_tensor_decomp(u,v,ones(n),ones(n)*rho);
% difuse
t = 100;
dt = .2;
niter=  round(t/dt);
kdisp = round(linspace(0,niter,5)); kdisp(1) = [];
f1 = rand(n);
v = []; k =1;
for i=1:niter
    f1 = f1 + dt * div( TensorMult(T, grad(f1,options) ) );
    if mod(i,50)==1
        f1 = perform_hist_eq(f1,'linear');
    end
    if i==kdisp(k)
        subplot(2,2,k);
        imageplot(perform_hist_eq(f1, 'linear'));
        k = k+1;
    end
end
