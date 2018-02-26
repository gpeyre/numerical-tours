W = matrix(0, nrow=p, ncol=k)
Elist = c()
tau = 0.01
niter = 1000
for (i in 1:niter)
{
    W = W - tau * nablaE(W)
    Elist = c(Elist, E(W))
}

plot(c(1: niter), Elist, type="l" ,col="blue", xlab="Iteration")