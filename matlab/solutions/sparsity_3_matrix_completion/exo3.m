niter = 200;
gamma = 1;
x = zeros(n);
E = [];
for i=1:niter
    x = ProxG( x - gamma * PhiS(Phi(x)-y), gamma*lambda );
    E(end+1) = 1/2*norm(Phi(x)-y, 'fro')^2 + lambda*sum(svd(x));
end
clf;
h = plot(E); axis tight; 
set(h, 'LineWidth', 2);
