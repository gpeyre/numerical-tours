




gD <- grad(phi, order=2)
d <- pmax(eps*array(1,c(n,n)), sqrt(apply(gD**2, c(1,2), sum)))
g <- gD/array(rep(d,2), c(dim(d),2))
G <- - W*d*div(g[,,1], g[,,2], order=2) - apply(gW*gD, c(1,2), sum)