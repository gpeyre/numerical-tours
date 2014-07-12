t = 40;
dt = .2;
niter=  round(t/dt);
kdisp = round(linspace(0,niter,5)); kdisp(1) = [];
rho_list = 1./[400 100 20 10];
for k=1:4
    rho = rho_list(k);
    T1 = perform_tensor_decomp(e1,e2,ones(n),ones(n)*rho);
    f1 = rand(n);
    for i=1:niter
        f1 = f1 + dt * div( TensorMult(T, grad(f1,options) ) );
        f1 = perform_hist_eq(f1,'linear');
    end
    imageplot(f1, ['1/\rho=' num2str(round(1/rho))], 2,2,k);
end
