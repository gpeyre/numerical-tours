niter = 200
Elist = rep(0, niter)
Clist = rep(0, niter)

for (i in 1:niter)
{
    # Update
    fold = f
    g = ProxFS(g + sigma * K(f1), sigma)
    f = ProxG(f - tau * KS(g), tau)
    f1 = f + theta * (f - fold)
    # Monitor the decay of energy
    Elist[i] = F(K(f))
    Clist[i] = snr(f0, f)
}
    
options(repr.plot.width=5, repr.plot.height=4)

plot(Elist, type="l", col=4, xlab="Iterations", ylab="E", main="Decay of the Energy")