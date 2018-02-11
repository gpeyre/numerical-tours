niter = 400
Elist = rep(0, niter)
Flist = rep(0, niter)

for (i in 1:niter)
{
    u = proxG(u - gamma * nablaF(u))
    x = y - As(u)
    Elist[i] = E(x)
    Flist[i] = F(u)
}
    
options(repr.plot.width=6, repr.plot.height=6)

plot(Elist[7:niter], type="l", ylim=c(450,600), col=4, xlab="", ylab="")
lines(-Flist[7:niter], col=3)

legend("topright", legend=c("E", "-F"), col=c(4,3), pch="-")