




J <- array(-(1/n), c(n,n))
diag(J) <- diag(J) + 1
K <- (-1/2.)* J %*% (D**2 %*% J)

ev <- eigen(K)
val <- ev$values ; Xstrain <- ev$vectors
ordering <- order(-Re(val))
val <- val[ordering] ; Xstrain <- Xstrain[,ordering]
val <- val[1:2] ; Xstrain <- Xstrain[,1:2]

Xstrain <- Xstrain * t(array(sqrt(val), c(2,n)))
Xstrain <- Re(t(Xstrain))

#plot points
plot(Xstrain[1,], Xstrain[2,], axes=FALSE, ann=FALSE, pch=19, col=color_function(X[1,], X[2,], X[3,]))

#plot vertices
I <- as.vector(B) ; J <- as.vector(t(NN))
xx <- array(0, c(2, k*n)) ; xx[1,] <- Xstrain[1,I] ; xx[2,] <- Xstrain[1,J]
yy <- array(0, c(2, k*n)) ; yy[1,] <- Xstrain[2,I] ; yy[2,] <- Xstrain[2,J]


for (i in 1:length(I)){
  lines(xx[,i], yy[,i], col="black", lw=1.5)
}


