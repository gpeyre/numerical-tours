h = zeros(n,1);
Flist = [];
tau = .5;
niter = 2000;
for i=1:niter
    h = h - tau * nablaF(h,K,y);
    Flist(i) = F(h,K,y);
end
clf;
plot(1:niter, Flist, 'LineWidth', 2);
axis tight;
SetAR(1/2); 
