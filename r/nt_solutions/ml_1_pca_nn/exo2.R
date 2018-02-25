library(grid)

options(repr.plot.width=5, repr.plot.height=5)

B = max(max(abs(Z[,1:2])))
q = 200
r = linspace(-B, B, q)
V1 = meshgrid(r)$X
U1 = meshgrid(r,r)$Y

z1 = cbind(c(U1), c(V1))
# test for different R
Rlist = c(1, 5, 10, 40)
col = c(1, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 0, 1, 1, 1, 0, 1, .5, .5, .5, 1, .5, .5, .5, 1)
col = matrix(col, nrow=3, ncol=(length(col) /3))

for (ir in 1:length(Rlist))
{
    R = Rlist[ir]
    D = distmat(Z[,1:2], z1)
    I = apply(D, 2, order)
    ys = y[I]
    dim(ys) = c(dim(Z)[1], dim(z1)[1])

    if (R==1)
    {
        C = ys[1,]
    }
    else
    {
        h = hist(ys[1:R,], 1:k)
        C = apply(h, 1, which.max)
    }

    # Need to plot two times in R
    for (i in 1:k)
    {
        
        I = (y==i)
        if (i == 1)
        {
            plot(Z[I,1], Z[I,2], xlim=c(-B, B), ylim=c(-B,B),
                 col=1 + i, pch=19, bty="n", xaxt="n", yaxt="n", xlab="", ylab="", main=paste("R=", R))
        }
        else
        {
           points(Z[I,1], Z[I,2], xlim=c(-B, B), ylim=c(-B,B), 
                  col=1 + i, pch=19, , bty="n", xaxt="n", yaxt="n",xlab="", ylab="", main=paste("R=", R)) 
        }
    }

    C = matrix(C, q, q)
    # maps class to color
    Cr = array(0, dim=c(q,q,3))
    for (i in 1:k)
    {
        for (a in 1:3)
        {
            Cr[,,a] = Cr[,,a] + (C==i) * col[a,i]
        }
    }
    # display
    grid.raster(aperm(Cr,c(2, 1, 3)), gp=gpar(alpha=0.5))
    
    for (i in 1:k)
    {
        I = (y==i)
        points(Z[I,1], Z[I,2], xlim=c(-B, B), ylim=c(-B,B), col=1 + i, 
               pch=19, xlab="", ylab="", bty="n", xaxt="n", yaxt="n", main=paste("R =", R))
    }
}