



r <- .95*n/2
p <- 128 # number of points on the curve
theta <- t( seq(0,2*pi,length=p+1) )
theta <- theta[1:length(theta)-1]
gamma0 <- n/2 * (1 + 1i) +  r*(cos(theta) + 1i*sin(theta))
gamma <- gamma0

imageplot(t(f))
cplot(gamma, type="l", pch=1, lty=1, lw=2, col="red", add=TRUE)