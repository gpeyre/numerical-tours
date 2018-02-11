nbiter = 500
x = y
E = c()
for (i in 1:nbiter)
{
    E = c(E, f(x, y, epsilon))
    x = x - tau * Gradf(x, y, epsilon)
}

plot(c(1:nbiter), E, xlab="Iterations", ylab="Energy (f)", main="Decay of the energy", col="blue", type="l")


options(repr.plot.width=5, repr.plot.height=5)
imageplot(clamp(x))