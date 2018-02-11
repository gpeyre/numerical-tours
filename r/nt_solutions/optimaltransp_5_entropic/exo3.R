



b <- rep(1, N)
niter <- 2000
Err_p <- c()
Err_q <- c()

for (i in 1:niter){
  a <- p/(xi%*%b)
  Err_q <- c(Err_q, norm(b*(xi%*%a) - q)/norm(q))
  b <- q /(t(xi)%*%a)
  Err_p <- c(Err_p, norm(a*(t(xi)%*%b) - p)/norm(p))
}

par(mfrow=c(2,1))

plot(log(Err_p), type="l", lwd = 2, col="blue", main="||pi-p||", xlab="", ylab="")

plot(log(Err_q), type="l", lwd = 2, col="blue", main="||pi^T-q||", xlab="", ylab="")