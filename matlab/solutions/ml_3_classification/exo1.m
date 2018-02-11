niter = 5000;
w = zeros(p+1,1);
Elist = [];
for i=1:niter
    w = w - tau * nablaE(w,AddBias(X),y);
    Elist(i) = E(w,AddBias(X),y);
end
ndisp = niter/2;
clf;
subplot(2,1,1);
plot(1:ndisp, Elist(1:ndisp), 'LineWidth', 2); axis tight;
% SetAR(1);
title('E(w_l)');
subplot(2,1,2);
plot(1:ndisp, log10(Elist(1:ndisp)-min(Elist)), 'LineWidth', 2); axis tight;
% SetAR(1);
title('log(E(w_l) - min E)');
