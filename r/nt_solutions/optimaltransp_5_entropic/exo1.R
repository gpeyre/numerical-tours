



b <- rep(1, N[2])
niter <- 300
Err_p <- c()
Err_q <- c()

for (i in 1:niter){
  a <- p/(xi%*%b)
  Err_q <- c(Err_q, norm(b*(t(xi)%*%a) - q)/norm(q))
  b <- q /(t(xi)%*%a)
  Err_p <- c(Err_p, norm(a*(xi%*%b) - p)/norm(p))
}

par(mfrow=c(2,1))

plot(log(Err_p + 1e-5), type="l", lwd = 2, col="blue", main="||pi-p||", xlab="", ylab="")

plot(log(Err_q + 1e-5), type="l", lwd = 2, col="blue", main="||pi^T-q||", xlab="", ylab="")