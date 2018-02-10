g = function(C,I){apply(C[,I], 1, sum)}

dlist = c(1:(N/20))
criter = matrix(0, length(dlist),3)

for (i in 1:length(dlist))
{
    s = twosparse(dlist[i])
    I = supp(s)
    criter[i,] = c(F(Phi, s), erc(Phi,I), werc(Phi,I))
}
    
criter[criter < 0] = NA

options(repr.plot.width=7, repr.plot.height=5)

matplot(dlist, criter, xlim=c(1, max(dlist)), ylim=c(min(criter), max(criter)), type="l", col=c(4, 3, 2),
        lty=1, xlab='', ylab='')
legend("right", legend=c('F', 'ERC', 'w-ERC'), 
       col=c(4, 3, 2), pch="-")
par(new=TRUE)
plot(dlist, dlist * 0 + 1, xlim=c(1, max(dlist)), ylim=c(min(criter), max(criter)), type="l", lty=2, ylab='', xlab='')
           