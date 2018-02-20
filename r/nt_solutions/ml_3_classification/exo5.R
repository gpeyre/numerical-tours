sigma_list = c(.1, .5, 1, 4)
niter = 4000
options(repr.plot.width=8, repr.plot.height=8)

par(mfrow=c(2,2))

for (is in 1:length(sigma_list))
{
    sigma = sigma_list[is]
    # grad descent
    K = kappa(X,X,sigma)
    Flist = c()
    tau = .5
    if (is == 4)
    {
        tau = .05
    }
    h = rep(0, n)

    for (i in 1:niter)
    {
        h = h - tau * nablaF(h,K,y)
        Flist = c(Flist, F(h,K,y))
    }

    # evaluate on a grid
    Theta = theta(kappa(G,X,sigma) %*% h)
    dim(Theta) = c(q, q)
    # Display the classification probability.
    image(t,t, Theta, xlab="", ylab="", col=color(10), xaxt="n", yaxt="n")
    par(new=TRUE)
    for (i in c(-1, 1))
    {
        I = (y==i)
        plot(X[I,1], X[I,2], col=(i + 3), xlim=c(min(X[,1]), max(X[,1])),  
             ylim=c(min(X[,2]), max(X[,2])), xlab="", ylab="", pch=16, main=paste("sigma =", sigma))
        par(new=TRUE)
    }
    par(new=FALSE)

}