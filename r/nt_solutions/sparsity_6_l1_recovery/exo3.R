N = 600
P = N/2
Phi = PhiRand(P, N)
klist = round(seq(1, P/7., length=20))
ntrials = 60
proba = matrix(0, length(klist), 3)

for (i in 1:length(klist))
{
    proba[i,] = 0
    k = klist[i]

    for (j in 1 : ntrials)
    {
        s = rep(0, N)
        I = sample(N)
        I = I[1:k]
        s[I] = sign(randn(n=k, m=1))
        proba[i, 1] = proba[i, 1] + (F(Phi, s) < 1)
        proba[i, 2] = proba[i, 2] + (erc(Phi, I) < 1)
        proba[i, 3] = proba[i, 3] + (werc(Phi, I) > 0) * (werc(Phi, I) < 1)
    }
}   

options(repr.plot.width=7, repr.plot.height=5)
matplot(klist, proba/ntrials, xlab="k", type="l", col=c(4, 3, 2), lty=1)
legend("right", legend=c('F <1', 'ERC <1', 'w-ERC <1'), 
       col=c(4, 3, 2), pch="-")