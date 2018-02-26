nbiter = 20
# Storing the results in a matrix
X = matrix(0, 2, nbiter + 1)
E = rep(0, nbiter + 1)

x = c(0.5,0.5)  # initial estimate of the solution
# Adding the first elements to the matrices
X[,1] = x
E[1] = f(x)

for (k in 1:nbiter){
    x = x - tau * GradF(x)
    X[,k + 1] = x 
    E[k + 1] = f(x)
    }


plot(c(1: (nbiter + 1)), log10(E), xlab="Iterations", ylab="log10(E)",
     main="Decay of the energy", col="blue", pch=15, type="l")

options(repr.plot.width=7, repr.plot.height=3.5)
filled.contour(t,t,t(F),nlevels=35, color.palette=topo.colors, main="Convergence of the gradient descent algorithm",
               plot.axes=lines(x=X[1,], y=X[2,], col="red", pch=15))