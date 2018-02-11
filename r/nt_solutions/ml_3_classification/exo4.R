options(repr.plot.width=5, repr.plot.height=4)

h = rep(0, n)
Flist = c()
tau = 0.5
niter = 2000

for (i in 1:niter)
{
    h = h - tau * nablaF(h,K,y)
    Flist = c(Flist, F(h,K,y))
}

plot(1:niter, Flist, type="l", col=4, xlab="", ylab="")