library (pracma)
library(akima)

perform_image_similitude <- function(M, mapping_type, u, u1, v, v1, w, w1){
  #
  #      perform_image_similitude
  
  #M1 = perform_image_similitude(M,mapping_type,u,u1,v,v1,w,w1);
  
  #Compute the affine similitude that map u to u1
  #and v to v1, and then resample the image M.
  #p and p1 are assumed to be in [0,1]=
  
  #If mapping_type=='similitude', compute a true similitude
  #T(x,y) = [a -b] * [x] + [c]
  #[b  a]   [y]   [d]
  #Solve the equations T(u)=u1 and T(v)=v1.
  
  #If mapping_type=='similitude', compute a true similitude
  #T(x,y) = [a  b] * [x] + [e]
  #[c  d]   [y]   [f]
  #Solve the equations T(u)=u1 and T(v)=v1 and T(w)=w1.
  
  #Copyright (c) 2006 Gabriel PeyrÈ
  
  if (mapping_type == "affine"){
    # the matrix of the linear system
    A <- t(Matrix(c(u[1], u[2], 0, 0, 1, 0,
           0, 0 , u[1], u[2], 0, 1,
           v[1], v[2], 0, 0, 1, 0,
           0, 0 , v[1], v[2], 0, 1,
           w[1], w[2], 0, 0, 1, 0,
           0, 0 , w[1], w[2], 0, 1), nrow=6, ncol=6))
    # the right hand size
    rhs <- c(u1,v1,w1)
    # solve
    z <- solve(A, rhs)
    # the similitude
    Q <- rbind(cbind(z[1],z[2]),cbind(z[3],z[4]))
    # the translation
    t <- c(z[5], z[6])
  } else {
    print ("Unknown mapping")
  }
  
  ### perform resampling ###
  
  # original grid in the warped domain
  n <- dim(M)[1]
  x <- seq(from=0, to=1, by=1/(n-1))
  X <- meshgrid(x,x)$X
  Y <- meshgrid(x,x)$Y
  # inverse warping P1=T^-1(P)
  P <- rbind(as.vector(X), as.vector(Y))
  P[1,] <- P[1,]-t[1]
  P[2,] <- P[2,] - t[2]
  P1 <- solve(Q) %*% P
  # reshape the results
  X1 <- P1[1,]
  Y1 <- P1[2,]
  
  M1 <- Matrix( bicubic(x, x, M, X1, Y1)$z, nrow=n, ncol=n)
  return (M1)
}