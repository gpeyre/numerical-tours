eta_list = seq(.1, 1, length=10)
ntrials = 20
mu_mean = c()
mu_std = c()
N = 500

for (i in 1:10)
{
    eta = eta_list[i]
    P = round(eta * N)
    c = c()
    
    for (k in 1:ntrials)
    {
        c = c(c, mu(PhiRand(P, N)))
    }    
    mu_mean = c(mu_mean, mean(c))
    mu_std = c(mu_std, std(c))
}

k_mean = 0.5 * (1 + 1./mu_mean)

options(repr.plot.width=5, repr.plot.height=3.5)
plot(eta_list, mu_mean, col=4, type="l", xlab="eta", ylab="mu")
plot(log10(eta_list), log10(k_mean), col=4, type="l", xlab="eta", ylab="std")