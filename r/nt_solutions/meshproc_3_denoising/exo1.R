klist <- c(1,2,4,8)
i <- 0
f1 <- f

for (k in (1:max(klist))){
  f1 <- tW %*% f1
  if (k == klist[i]){
    my_cmap <- cbind(f1,f1,f1)
    scatterplot3d(X0[1,], X0[2,], X0[3,], axis=FALSE, grid=FALSE, box=FALSE, type="p", pch=19, color=my_cmap[,1], scale.y=0)
    i <- i+1 
    }
}
