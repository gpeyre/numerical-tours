options(repr.plot.width=5, repr.plot.height=3)
D = distmat(X0, X1)
I = apply(D, 2, order)
ys = y0[I]
dim(ys) = c(dim(X0)[1], dim(X1)[1])
Rmax = 50
S = c()

for (R in 1:Rmax)
{
    if (R==1)
    {
        C = ys[1,]
    }
    else
    {
        h = hist(ys[1:R,], 1:k)
        C = apply(h, 1, which.max)
    }
    # Correct classfication
    S = c(S, sum(C==y1) / n1)      
}

barplot(S, col=4, ylim=c(min(S), 1), names.arg=c(1:Rmax), xlab="R",xpd=FALSE)