s = 31
sel = sample(c(1 : n))
x0 = zeros(n, 1)
x0[sel[1:s]] = 1
y = A %*% x0

tx = zeros(n,1)

for (i in 1:niter)
{
    tx = (1 - (mu/2)) * tx + (mu / 2) * rproxG(rproxF(tx, y), gamma)
    x = proxF(tx,y)
}


options(repr.plot.width=5, repr.plot.height=3)
stemplot(x0, ylab="", xlab="", main="Original Signal", ylim=c(-1,1))
stemplot(x, ylab="", xlab="", main="Recovered Signal", ylim=c(-1,1))