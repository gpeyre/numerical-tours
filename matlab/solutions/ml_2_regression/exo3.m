Jlist = [];
niter = 400;
w = zeros(p,1);
for i=1:niter
    Jlist(end+1) = J(w,lambda);
    w = ISTA(w,lambda,tau);
end
ndisp = niter/4;
clf;
subplot(2,1,1);
plot(1:ndisp,Jlist(1:ndisp), 'LineWidth', 2);
axis tight; 
title('J(w_k)');
subplot(2,1,2);
e = log10(Jlist(1:ndisp)-min(Jlist));
plot(1:ndisp,e-e(1), 'LineWidth', 2); 
axis tight;
title('log(J(w_k)-min J)');
