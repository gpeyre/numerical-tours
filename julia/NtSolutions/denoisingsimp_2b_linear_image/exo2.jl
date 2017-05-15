mulist = linspace(.1, 3.5, 31)
# err = arrayfun(@(mu)norm(x0-denoise(y,mu), 'fro' ), mulist);
err = zeros(length(mulist))
for i in 1 : length(mulist)
	mu = mulist[i]
	err[i] = norm(x0 - denoise(y, mu))
end
clf
#h1, = plot(mulist,err); axis("tight");
plot(mulist,err); axis("tight")
# set_label('\mu', '|y-x_0|');
# retrieve the best denoising result
i = mapslices(indmin, err, 1)[1]
mu = mulist[i]
