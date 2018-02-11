mulist = linspace(.1,3.5,31);
err = arrayfun(@(mu)norm(x0-denoise(y,mu), 'fro' ), mulist);
clf;
h = plot(mulist,err); axis('tight');
set_label('\mu', '|y-x_0|');
% retrieve the best denoising result
[~,i] = min(err);
mu = mulist(i);
