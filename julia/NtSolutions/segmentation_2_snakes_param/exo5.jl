G = grad(ff)
G = sqrt(sum(G.^2, 3))
sigma = 1.5
G = gaussian_blur(G[2 : end, 2 : end, 1], sigma)
G = min(G, .4)
W = rescale(-G, .4, 1)
clf
imageplot(W)