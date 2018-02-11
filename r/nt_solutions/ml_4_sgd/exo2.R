niter = 100000
nsamples = 10
err_rate = 50
ElistS = matrix(0, nrow=(niter / err_rate), ncol=nsamples)

for (is in 1:nsamples)
{
    w = rep(0, p)
    for (it in 1:niter)
    {
        if ((it %% err_rate) == 1)
        {
            ElistS[1 + floor((it-1)/err_rate), is] = E(w,X,y)
        }
        tau = tau0 / (1 + (it/l0))
        i = (1 + floor(rand() * n))[1] # draw uniformly
        w = w - tau * nablaEi(w,i)
    }
}

# First plot
matplot(seq(1, niter, err_rate), ElistS, type="l", col=4, lwd=2, xlab="", ylab="", main="E(wi)")
lines(1 + c(0 : (length(Elist) - 1)) * n, Elist, type="l", lty=2)

u = log10(ElistS - min(Elist))
v = log10(Elist - min(Elist))
matplot(seq(1, niter, err_rate), u, type="l", col=4, lwd=2, xlab="", ylab="", main="log(E(w_l) - min E)")
lines(1 + c(0 : (length(Elist) - 1)) * n, v, type="l", lty=2)