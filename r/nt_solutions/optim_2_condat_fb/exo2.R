nbiter = 4000

# Storing the results in a matrix
E = rep(0, nbiter + 1)
x = y
E[1] = 0.5 * norm(Phi(x) - y)^2 + Lambda * l1_norm(x)

for (k in 1:nbiter)
{
    x = proxg(x - gamma * grad_f(x), gamma, Lambda)
    E[k + 1] = 0.5 * norm(Phi(x) - y)^2 + Lambda * l1_norm(x)
}

E_plot = (E[1:(nbiter / 4)] - E[nbiter + 1]) / (E[1] - E[nbiter + 1])
plot(log10(E_plot), type="l", col="blue", ylab="log10((E - E*) / (E0 - E*))", xlab="i")

par(mfrow=c(1,2))
options(repr.plot.width=7, repr.plot.height=3.5)
stemplot(x0, col="blue", ylab="", xlab="", main="Signal x0")
stemplot(x, col="blue", ylab="", xlab="", main="Recovered Signal")