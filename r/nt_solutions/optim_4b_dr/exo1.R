lun = c()
err = c()
tx = zeros(n,1)

for (i in 1:niter)
{
    tx = (1 - (mu/2)) * tx + (mu / 2) * rproxG(rproxF(tx, y), gamma)
    x = proxF(tx,y)
    lun = c(lun, l1_norm(x))
    err = c(norm(y - A %*% x))
}

plot(lun, type="l", col="blue", ylab="||x|_1", xlab="Nb of iterations")

plot(log10(lun[1: (niter/2)]- lun[niter]), col="blue", ylab="log10(||x|_1)", xlab="Nb of iterations", type="l")


options(repr.plot.width=5, repr.plot.height=3)
stemplot(x0, ylab="", xlab="", main="Original Signal", ylim=c(-1,1))
stemplot(x, ylab="", xlab="", main="Recovered Signal", ylim=c(-1,1))