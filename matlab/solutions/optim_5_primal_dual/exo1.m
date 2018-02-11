niter = 200;
E = []; C = [];
for i=1:niter    
    % update
    fold = f;
    g = ProxFS( g+sigma*K(f1), sigma);
    f = ProxG(  f-tau*KS(g), tau);
    f1 = f + theta * (f-fold);
    % monitor the decay of the energy
    E(i) = F(K(f));
    C(i) = snr(f0,f);
end
clf;
h = plot(E);
set(h, 'LineWidth', 2);
axis('tight');
