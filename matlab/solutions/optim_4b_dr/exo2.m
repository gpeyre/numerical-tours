s = 31;
sel = randperm(n);
x0 = zeros(n,1); x0(sel(1:s))=1;
y = A*x0;
tx = zeros(n,1);
for i=1:niter
    tx = (1-mu/2)*tx + mu/2*rproxG( rproxF(tx,y),gamma );
end
x = proxF(tx,y);
clf;
subplot(2,1,1);
plot_sparse_diracs(x0);
set_graphic_sizes([], 15);
title('Original Signal');
subplot(2,1,2);
plot_sparse_diracs(x);
set_graphic_sizes([], 15);
title('Recovered by L1 minimization');
