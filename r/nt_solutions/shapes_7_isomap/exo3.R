



niter <- 150
stress <- c()
Xstress <- X
ndisp <- c(1, 5, 10, min(niter,100), Inf)
plt <- 1


par(mfrow=c(2,2))


for (i in 1:niter){

  if (ndisp[plt]==i){
    
    #swiss roll
    s3d <- scatterplot3d(Xstress[1,], Xstress[2,], Xstress[3,], axis=F, grid=F, box=F, type="p", pch=19, color=color_function(X[1,], X[2,], X[3,]))
    #graph
    xx <- array(0, c(2, k*n)) ; xx[1,] <- Xstress[1,I] ; xx[2,] <- Xstress[1,J]
    yy <- array(0, c(2, k*n)) ; yy[1,] <- Xstress[2,I] ; yy[2,] <- Xstress[2,J]
    zz <- array(0, c(2, k*n)) ; zz[1,] <- Xstress[3,I] ; zz[2,] <- Xstress[3,J]
    
    for (i in 1:length(I)){
      s3d$points3d(xx[,i], yy[,i], zz[,i], type="l", lw=1.5)
    }
    plt <- plt+1
  }
  
  # Compute the distance matrix.
  D1 <- array( rep(apply(Xstress**2, 2, sum), n), c(n,n) )
  D1 <- sqrt(pmax(D1 + t(D1) - 2*t(Xstress)%*%Xstress, 0))
  
  
  # Compute the scaling matrix.
  B <- - D / pmax(D1, array(1e-10, dim(D1)))
  diag(B) <- diag(B) - apply(B, 2, sum)
  
  # update
  Xstress <- t( B %*% t(Xstress) )/n
  
  
  # record stress
  stress <- c(stress, sqrt(sum(abs(D-D1)**2)/n**2))
  
}