


T <- 3*sigma
w <- 4
Mspin <- array(0, c(n,n,n))

for (i in 1:w**3){
  # shift the image
  MnoisyC <- circshift(Mnoisy, c(dX[i],dY[i],dZ[i]))
  # denoise 
  MW <- perform_haar_transf(MnoisyC, 1, + 1)
  MWT <- perform_thresholding(MW, T, "hard")
  M1 <- perform_haar_transf(MWT, 1, -1)
  # back
  M1 <- circshift(M1, c(-dX[i],-dY[i],-dZ[i]))
  # average the result
  Mspin <- Mspin*(i-1)/(i) + M1/(i)
}