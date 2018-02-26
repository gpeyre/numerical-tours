tau0 = 0.05
ell0 = 100
nsamples = 10
err_rate = 50
ElistA = matrix(0, nrow=(niter / err_rate), ncol=nsamples)

options(repr.plot.width=7, repr.plot.height=5)

for (is in 1:nsamples)
{
    w = rep(0, p)
    w1 = w
    for (it in 1:niter)
    {
        if ((it %% err_rate) == 1)
        {
            ElistA[1 + floor((it-1)/err_rate), is] = E(w,X,y)
        }
        tau = tau0 / (1 + sqrt((it-1)/ell0))
        i = (1 + floor(rand() * n))[1] # draw uniformly
        w1 = w1 - tau * nablaEi(w1,i)
        w = (1/it) * w1 + (1 - (1/it)) * w
    }
}

vmin = min(min(ElistS), min(ElistA), min(Elist))

# Getting the boundaries for the plot
min_1 = min(log10(ElistS-vmin))
min_2 = min(log10(ElistA-vmin))
max_1 = max(log10(ElistS-vmin))
max_2 = max(log10(ElistA-vmin))
min_plot = min(min_1, min_2)
max_plot = max(max_1, max_2)

matplot(seq(1, niter, err_rate), log10(ElistS-vmin+1e-18), type="l", col=4, lwd=2, 
        xlab="", ylab="", main="log(E(w_l) - min E)", ylim=c(min_plot, max_plot))
par(new=TRUE)
matplot(seq(1, niter, err_rate), log10(ElistA-vmin + 1e-18), type="l", col=2,
        lwd=2, xlab="", ylab="", ylim=c(min_plot, max_plot))
legend("topright", legend=c("SGD", "SGA"), col=c(4,2), pch="-")