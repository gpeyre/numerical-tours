tau = .5;
niter = 5000;
w = zeros(p,1);
Elist = [];
for i=1:niter
    w = w - tau * nablaE(w,X,y);
    Elist(i) = E(w,X,y);
end
ndisp = niter/4;
clf;
subplot(2,1,1);
plot(1:ndisp, Elist(1:ndisp), 'LineWidth', 2); axis tight;
title('E(w_l)');  set(gca, 'FontSize', fs);
subplot(2,1,2);
plot(1:ndisp, log10(Elist(1:ndisp)-min(Elist)), 'LineWidth', 2); axis tight;
title('log(E(w_l) - min E)'); set(gca, 'FontSize', fs);
