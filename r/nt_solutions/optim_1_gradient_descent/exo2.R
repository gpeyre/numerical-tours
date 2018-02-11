taulist = c(.3, 1, 1.7) / eta
xinit = list(c(.7, .7), c(-.7, .5), c(-.7, -.6))
nbiter = 100

listmat = rep(list(matrix(0, 2, nbiter + 1)), 3)

for (k in 1:length(taulist))
{
    tau = taulist[k]

    x = xinit[[k]]
    listmat[[k]][, 1] = x

    for (i in 1:nbiter){
        x = x - tau * GradF(x)
        listmat[[k]][,i + 1] = x
        }
}

options(repr.plot.width=7, repr.plot.height=3.5)
filled.contour(t,t,t(F),nlevels=35, color.palette=topo.colors, 
               main="Convergence according to the learning rate",
               plot.axes={lines(x=listmat[[1]][1,], y=listmat[[1]][2,], col="red", pch=15)
                          lines(x=listmat[[2]][1,], y=listmat[[2]][2,], col="green", pch=15)
                          lines(x=listmat[[3]][1,], y=listmat[[3]][2,], col="black", pch=15)})

legend(0,-10, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"))