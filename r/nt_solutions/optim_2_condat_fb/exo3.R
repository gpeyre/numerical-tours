nbiter = 4000

mu_list = c(-0.5, 0, 0.5, 0.8, 0.9, 0.95, 0.99)

# Storing the results
Erelax = rep(list(rep(0, nbiter + 1)), length(mu_list))

ind = 1

for (mu in mu_list)
{
    x = y
    z = y
    Erelax[[ind]][1] = 0.5 * norm(Phi(x) - y)^2 + Lambda * l1_norm(x)
    for (k in 1:nbiter)
    {
        xold = x
        zold = z
        x = proxg(z - gamma * grad_f(z), gamma, Lambda)
        z = x + mu * (x - xold)
        Erelax[[ind]][k + 1] = 0.5 * norm(Phi(x) - y)^2 + Lambda * l1_norm(x)
    }
    ind = ind + 1
}

options(repr.plot.width=8, repr.plot.height=5)

E_plot = (Erelax[[1]][1 : (nbiter / 4)] - Erelax[[1]][nbiter + 1]) / (Erelax[[1]][1] - Erelax[[1]][nbiter + 1]) 
plot(log10(E_plot), type="l", col=1, ylab="log10((E - E*) / (E0 - E*))", xlab="i", ylim=c(-20, 0))

for (i in 2:length(mu_list))
{
    E_plot = (Erelax[[i]][1 : (nbiter / 4)] - Erelax[[i]][nbiter + 1]) / (Erelax[[i]][1] - Erelax[[i]][nbiter + 1]) 
    lines(log10(E_plot), col=i)
}

vec_legend = lapply(mu_list, function(x){paste("mu = ", x)})

legend("topright", legend=vec_legend, col=1:length(mu_list), pch="-")