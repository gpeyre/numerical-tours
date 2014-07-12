T = 20;
% step size
tau = .2;
% number of iteration
niter = round(T/tau);
Mheat = perform_hist_eq(randn(n), x);
clf; k = 0;
for i=1:niter
    % compute the gradient
    G = div(grad(Mheat));
    % descent 
    Mheat = perform_hist_eq(Mheat + tau*G, x);
    if mod(i,floor(niter/4))==0
        k = k+1;
        imageplot(Mheat, strcat(['T=' num2str(T*k/4,3)]), 2,2,k );
    end
end
