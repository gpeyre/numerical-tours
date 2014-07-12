p = round(rho*n^2);
sel = randperm(n^2); sel = sel(1:p);
M = M0; 
M(sel) = M(sel) + randn(size(sel)) * sigma;
clf;
imageplot(M0, 'Original', 1,2,1);
imageplot(clamp(M), 'Noisy', 1,2,2);
