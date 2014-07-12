G = []; 
F = [];
x = x0;
x = zeros(n); 
tx = x;
niter = 500;
for i=1:niter
    tx = (1-mu/2)*tx + mu/2*rProxG( rProxF(tx,gamma),gamma );
    x = ProxF(tx,gamma);
    G(i) = sum(svd(x));
    F(i) = norm(y-Phi(x));
end
clf;
h = plot(G);
set(h, 'LineWidth', 2);
axis tight;
