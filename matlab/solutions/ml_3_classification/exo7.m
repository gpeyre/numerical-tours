W = zeros(p,k);
Elist = [];
tau = .01;
niter = 500;
for i=1:niter
    W = W - tau * nablaE(W);
    Elist(i) = E(W);
end
clf;
plot(1:niter, Elist, 'LineWidth', 2);
axis tight;
SetAR(1/2); 
