a = 10
b = 10
tx = seq(-a, a, length=q)
ty = seq(-b, b, length=q)
B = as.vector(meshgrid(ty, tx)$X)
A = as.vector(meshgrid(ty, tx)$Y)
G = matrix(c(A, B), nrow=length(A), ncol=2)
offs = c(0.3, 1, 3, 5)
niter = 10000

par(mfrow=c(2,2))
options(repr.plot.width=6, repr.plot.height=6)

for (io in 1:length(offs))
{
    # generate data
    omega = offs[io] * 1.5
    X = rbind(randn(n/2,2), randn(n/2,2) + rep(1, n/2) * omega)
    y = c(rep(1, n/2), rep(-1, n/2))
    # run gradient descent
    w = rep(0, p + 1)
    for (i in 1:niter)
    {
    w = w - tau * nablaE(w, AddBias(X), y)
    }

    # display
    Theta = theta(AddBias(G) %*% w)
    dim(Theta) = c(q, q)

    image(tx,ty, Theta, xlab="", ylab="", col=color(5), xaxt="n", yaxt="n")
    
    par(new=TRUE)
    for (i in c(-1, 1))
    {
        I = (y==i)
        plot(X[I,1], X[I,2], col=(i + 3), xlim=c(min(X[,1]), max(X[,1])), 
           ylim=c(min(X[,2]), max(X[,2])), xlab="", ylab="", pch=16, xaxt="n", yaxt="n")  
        par(new=TRUE)
    }
    par(new=FALSE)
    
}