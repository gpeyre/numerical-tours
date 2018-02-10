N = 600
P = N/2
Phi = PhiRand(P, N)
klist = round(seq(1, P/15., length=10))
ntrials = 2000
rip_val = matrix(0, length(klist), 2)

for (i in 1:length(klist))
{
    rip_val[i,1:2] = 0
    k = klist[i]
    for (j in 1:ntrials)
    {
        I = sample(N)
        I = I[1:k]
        tmp = ric(Phi[,I])
        a = tmp[1]
        b = tmp[2]
        rip_val[i,1:2] = pmax(rip_val[i,1:2] , c(a,b))    
    }
}

options(repr.plot.width=7, repr.plot.height=5)

matplot(klist, rip_val, type='l', lty=1, main=paste("N=", N, ", P=", P), xlab="",
        ylab="", ylim=c(min(rip_val), max(rip_val)), col=c(4, 3))
legend("right", legend=c('delta^2_k', 'delta^2_k', '0.41'), 
       col=c(4, 3, 2), pch="-")
par(new=TRUE)
plot(klist, klist*0 + sqrt(2)-1, xlab="", ylab="", type='l', lty=2, ylim=c(min(rip_val), max(rip_val)))