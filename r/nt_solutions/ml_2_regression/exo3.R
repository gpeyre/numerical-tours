niter = 400
w = rep(0, p)
Jlist = rep(0, niter)

for (i in 1:niter)
{
    Jlist[i] = J(w,Lambda)
    w = ISTA(w, Lambda,tau)
}

ndisp = niter/4
plot(Jlist[1: ndisp], main='J(w_k)',type="l", col=4, xlab="", ylab="")
e = log10(Jlist[1:ndisp] - min(Jlist) + 1e-20)
plot(e - e[1], type="l", col=4, main="log(J(w_k)-min J)", xlab="", ylab="")