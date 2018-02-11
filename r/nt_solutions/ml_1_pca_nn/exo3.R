options(repr.plot.width=6, repr.plot.height=6)

I = sample(n)
I = I[1:k]
C = X[I,]
niter = 16
it_dist = 0

# Variables for plotting
l_position_plots = list(c(0.1,0.4,0.6,0.9), c(0.1,0.4,0.1,0.4), c(0.6,0.9,0.6,0.9), c(0.6,0.9,0.1,0.4))
plot_iter = 1

par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))

for (it in 1:niter)
{
    # NN
    D = distmat(X,C)
    yb = apply(D, 1, which.min)
    # display
    if (it <= 3 || it == niter)
    {
        it_dist = it_dist + 1
        par(new=FALSE)
        if (plot_iter ==1)
        {
            par(fig=l_position_plots[[plot_iter]])
        }
        else
        {
            par(fig=l_position_plots[[plot_iter]], new=TRUE)
        }
        plot_iter = plot_iter + 1

        for (i in 1:k)
        {
            I = (yb == i)
            plot(Z[I,1], Z[I,2], col=(i+1), xlim=c(min(Z[,1]), max(Z[,1])),  
             ylim=c(min(Z[,2]), max(Z[,2])), xaxt='n', yaxt='n', 
             ylab="", xlab="", main=paste("Iter #", it, end=""), bty="n")
            par(new=TRUE)
            
        }
    CV = (C - rep(colMeans(X), rep.int(nrow(C), ncol(X)))) %*% V
    for (i in 1:k)
    {
        points(CV[i, 1], CV[i, 2], , xlim=c(min(Z[,1]), max(Z[,1])),  
             ylim=c(min(Z[,2]), max(Z[,2])), xaxt='n', yaxt='n', 
             ylab="", xlab="", lwd=5, bty="n")
    }
    }
    # update centroids
    for (l in 1:k)
    {
        C[l,] = apply(X[yb==l,], 2, mean)
    }
}


options(repr.plot.width=5, repr.plot.height=3)

for (l in (1:k))
{
    I = (yb==l)
    h = hist(y[I],1:k)
    h = h/sum(h)
    barplot(h, axes=FALSE, col="blue")

    # add y-axis
    axis(side = 2, pos = -0.2)
    # add x-axis with offset positions, with ticks, but without labels.
    axis(side = 1, at = 0:(k - 1), labels = FALSE)
    # add x-axis with centered position, with labels, but without ticks.
    axis(side = 1, at = 1:k - 0.5, tick = FALSE, labels = 1:k)
}