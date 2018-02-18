q = 20  # Number of paths computed in parallel.
niter = 1000

x = matrix(0, q, niter)
x[,0] = randn(q,1) - .5 # Initial conditions.

for (i in 2:niter)
{
    u = randn(q,1) > 0.5
    xx = x[, i - 1]
    tau = 1.0 / (10 + i)
    x1 = xx - tau * (u * df(xx,1) + (1-u) * df(xx,2))
    x[,i] = x1
}


matplot(t(x), col='red', type="l", xlab="", ylab="")
