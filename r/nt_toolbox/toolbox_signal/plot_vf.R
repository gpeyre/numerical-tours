



plot_vf <- function(velocities){
  ####
  # velocities is supposed to be of shape nxnx2
  ####
  n_v <- dim(velocities)[1]
  u <- velocities[,,1]
  v <- velocities[,,2]
  grid <- meshgrid_2d(0:(n_v-1), 0:(n_v-1))
  x <- grid$X ; y <- grid$Y
  plot(0,0,xlim=c(0,n_v), ylim=c(0,n_v), type="n", axes=FALSE, ann=FALSE)
  quiver(x,y,u,v, scale=0.7,length=0.07, angle=15, col="blue", lwd=1, code=1)
  
}


