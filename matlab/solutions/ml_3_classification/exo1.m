w = zeros(n,1);
Elist = [];
tau = .5;
niter = 500;
for i=1:niter
    w = w - tau * nablaE(w);
    Elist(i) = E(w);
end
ndisp = 200;
clf;
subplot(2,1,1);
plot(1:ndisp, Elist(1:ndisp), 'LineWidth', 2); axis tight;
% SetAR(1); 
title('E(w_l)');
subplot(2,1,2);
plot(1:ndisp, log10(Elist(1:ndisp)-min(Elist)), 'LineWidth', 2); axis tight;
% SetAR(1); 
title('log(E(w_l) - min E)');
