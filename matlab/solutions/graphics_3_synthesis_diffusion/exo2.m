T = 3;
% step size
tau = .005;
% avoid instabilities
epsilon = 1e-5;
% number of iteration
niter = round(T/tau);
Mtv = perform_hist_eq(randn(n), x);
clf; k = 0;
for i=1:niter
    % compute the gradient
    G = grad(Mtv);
    d = max(epsilon, sqrt(sum(G.^2,3)));
    G = div( G ./ repmat(d, [1 1 2]) );
    % descent 
    Mtv = perform_hist_eq( Mtv + tau*G, M);
    if mod(i,floor(niter/4))==0
        k = k+1;
        Gr = grad(Mtv);
        tv = sum(sum( sqrt(sum(Gr.^2,3)), 2 ), 1);
        imageplot(Mtv, strcat(['T=' num2str(T*k/4,3) ', TV=' num2str(tv)]), 2,2,k );
    end
end
