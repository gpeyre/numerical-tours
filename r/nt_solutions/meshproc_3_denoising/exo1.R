klist <- c(1,2,4,8)
i <- 1
f1 <- f

for (k in (1:max(klist))){
  f1 <- tW %*% f1
  if (k == klist[i]){
    my_cmap <- cbind(f1,f1,f1)
    points3D(matrix(X0[1,]), matrix(X0[2,]), matrix(X0[3,]),cex=0.05, axis=FALSE, grid=FALSE, box=FALSE,colkey=FALSE, type="p", pch=20, col=as.color(f1))
    i <- i+1 
    }
}
