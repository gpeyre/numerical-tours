c = [];
q = 100;
niter = 60;
for i=1:niter
    c = [c C(repmat(x0,[1 q]),randn(P,q))];
end
sigma0 = std(c);
t = linspace(-4*sigma0,4*sigma0,31);
h = hist(c, t);
clf;
bar(t,h); axis('tight');
