tx = x0
niter = 2000
for (i in 1:niter)
{
    tx = (1 - (mu / 2)) * tx + (mu / 2) * rproxG(rproxF(tx,y),gamma)
    x = proxF(tx,y)
}

proba = c()
E = 1 * (apply(abs(x - x0),2,mean) < .05)

for (j in c(1:length(slist)))
{
    s = slist[j]
    proba = c(proba, mean(E[Slist==s]))
}

plot(slist, proba, type="l", col="blue", xlab="Sparsity",
     ylab="Proba", main="Evaluation of the CS Recovery Probability")