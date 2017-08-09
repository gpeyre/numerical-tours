h = zeros(p,1);
Flist = [];
tau = 1000;
niter = 200;
for i=1:niter
    h = h - tau * nablaF(h);
    Flist(i) = F(h);
end
clf;
plot(1:niter, Flist, 'LineWidth', 2);
axis tight;
SetAR(1/2); 
