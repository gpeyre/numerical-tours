



G <- grad(f)
G <- sqrt(apply(G**2, c(1,2), sum))
sigma <- 1.5
G <- gaussian_blur(G,sigma)
G <- pmin(G,.4)
W <- rescale(-G,.4,1)

imageplot(W)