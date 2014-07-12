Tlist = linspace(.8,1.2,15)*sigma;
MW = perform_curvelet_transform(M, options);
err = [];
for i=1:length(Tlist)
    MWT = perform_thresholding(MW, Tlist(i), 'hard');
    M1 = perform_curvelet_transform(MWT, options);
    err(end+1) = snr(M0,M1);
end
clf;
h = plot(Tlist/sigma,err);
set(h, 'LineWidth', 2);
axis('tight');
[tmp,i] = max(err);
MWT = perform_thresholding(MW, Tlist(i), 'hard');
Mcurv = perform_curvelet_transform(MWT, options);
