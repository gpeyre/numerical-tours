gamma = 1 / L
nbiter = 4000
x = y
z= y
t = 1

E_fista = rep(0, nbiter + 1)
E_fista[1] = 0.5 * norm(Phi(x) - y)^2 + Lambda * l1_norm(x)

for (k in 1 : nbiter)
{
    xold = x
    told = t
    x = proxg(z - gamma * grad_f(z), gamma, Lambda)
    t = (1 + sqrt(1 + 4 * t^2) ) / 2;
    z = x + (told - 1) / t * ( x - xold );
    E_fista[k + 1] = 0.5 * norm(Phi(x) - y)^2 + Lambda * l1_norm(x)  
}

# Taking the minimum value obtained so far
E_min = min(min(E_fista), min(Erelax[[7]]), min(E))


# Vector to plot
E_fista_plot = (E_fista[1 : (nbiter / 4)] - E_min) / (E_fista[1] - E_min)
E_relax_plot = (Erelax[[7]][1 : (nbiter / 4)] - E_min) / (Erelax[[7]][1] - E_min)
E_plot = (E[1 : (nbiter / 4)] - E[nbiter]) / (E[1] - E[nbiter])

options(repr.plot.width=7, repr.plot.height=6)

plot(log10(E_fista_plot), type="l", col=1, ylim=c(-12,0), ylab="log10((E - E*) / (E0 - E*))")
lines(log10(E_relax_plot), col=2)
lines(log10(E_plot), col=3)

legend("top", legend=c("Fista", "FB-relax, mu=0.99", "FB"), 
       col=1:3, pch="-")