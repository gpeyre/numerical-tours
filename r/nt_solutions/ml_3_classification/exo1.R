options(repr.plot.width=5, repr.plot.height=4)
niter = 5000
w = rep(0, p + 1)
Elist = c()

for (i in 1:niter)
{
    w = w - tau * nablaE(w, AddBias(X), y)
    Elist = c(Elist, E(w,AddBias(X),y))
}

ndisp = niter/2
plot(1:ndisp, Elist[1:ndisp], type='l', col=4, xlab="", ylab="",  main="E(w_i)")
plot(1:ndisp, log(Elist[1:ndisp]-min(Elist), 10), type='l', col=4,
     xlab="", ylab="",  main="log(E(w_l) - min E)")

