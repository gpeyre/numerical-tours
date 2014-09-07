mulist = linspace(.1,3.5,31);
# err = arrayfun(@(mu)norm(x0-denoise(y,mu), 'fro' ), mulist);
err = zeros([len(mulist),1]) 
for i in arange(0,len(mulist)):
	mu = mulist[i];
	err[i] = norm(x0-denoise(y,mu));
clf;
h1, = plot(mulist,err); axis('tight');
# set_label('\mu', '|y-x_0|');
# retrieve the best denoising result
i = argmin(err);
mu = mulist[i];
