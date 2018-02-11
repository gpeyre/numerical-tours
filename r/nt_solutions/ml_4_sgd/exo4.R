tau = .002 / n
nsamples = 10
err_rate = 50

ElistG = matrix(0, nrow=(niter / err_rate), ncol=nsamples)

options(repr.plot.width=7, repr.plot.height=5)

for (is in 1:nsamples)
{
    w = rep(0, p)
    G = matrix(0, nrow=p, ncol=n)
    g = rep(0, p)
    for (it in 1:niter)
    {
        if ((it %% err_rate) == 1)
        {
            ElistG[1 + floor((it-1)/err_rate), is] = E(w,X,y)
        }
        i = (1 + floor(rand() * n))[1] # draw uniformly
        g1 = nablaEi(w,i)
        #update grad
        g = g - G[, i] + g1
        G[,i] = g1
        w = w - tau * g
    }
}


# Getting the boundaries for the plot
min_1 = min(log10(ElistS-min(Elist)))
min_2 = min(log10(ElistA-min(Elist)))
min_3 = min(log10(ElistG-min(Elist)))
max_1 = max(log10(ElistS-min(Elist)))
max_2 = max(log10(ElistA-min(Elist)))
max_3 = max(log10(ElistG-min(Elist)))
min_plot = min(min_1, min_2, min_3)
max_plot = max(max_1, max_2, max_3)

matplot(seq(1, niter, err_rate), log10(ElistS-min(Elist)), type="l", col=4, lwd=2, 
        xlab="", ylab="", main="log(E(w_l) - min E)", ylim=c(min_plot, max_plot))
par(new=TRUE)
matplot(seq(1, niter, err_rate), log10(ElistA-min(Elist)), type="l", col=2,
        lwd=2, xlab="", ylab="", ylim=c(min_plot, max_plot))
par(new=TRUE)
matplot(seq(1, niter, err_rate), log10(ElistG-min(Elist)), type="l", col=3,
        lwd=2, xlab="", ylab="", ylim=c(min_plot, max_plot))
legend("topright", legend=c("SGD", "SGA", "SAG"), col=c(4,2,3), pch="-")