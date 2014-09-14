G = grad(f)
G = sqrt(sum(G**2,2))
sigma = 1.5
G = gaussian_blur(G,sigma)
G = minimum(G,.4)
W = rescale(-G,.4,1)
clf
imageplot(W)
