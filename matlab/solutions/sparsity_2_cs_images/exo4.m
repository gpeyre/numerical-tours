G = []; 
F = [];
tx = zeros(N,1);
niter = 300;
for i=1:niter
    tx = (1-mu/2)*tx + mu/2*rProxG( rProxF(tx,gamma),gamma );
    x = ProxF(tx,gamma);
    G(i) = norm(x,1);
    F(i) = norm(y-Phi(x));
end
clf;
h = plot(G);
set(h, 'LineWidth', 2);
axis tight;
