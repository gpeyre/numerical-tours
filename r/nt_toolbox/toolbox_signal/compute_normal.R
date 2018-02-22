crossp <- function(x,y){
  # x and y are (m,3) dimensional
  z <- x
  z[1,] = x[2,]*y[3,] - x[3,]*y[2,]
  z[2,] = x[3,]*y[1,] - x[1,]*y[3,]
  z[3,] = x[1,]*y[2,] - x[2,]*y[1,]
  return (z)
}

compute_normal <- function(X0, F){
  
  #compute_normal - compute the normal of a triangulation
  
  #[normal,normalf] = compute_normal(vertex,face);
  
  #normal(i,:) is the normal at vertex i.
  #normalf(j,:) is the normal at face j.
  eps <- .Machine$double.eps
  nface <- dim(F)[2]
  nvert <- dim(X0)[2]
  normal <- matrix(0, nrow=3, ncol = nvert)
  
  # unit normals to the faces
  normalf <- crossp(X0[,(F+1)[2,]] - X0[,(F+1)[1,]],X0[,(F+1)[3,]] - X0[,(F+1)[1,]])
  d <- sqrt(colSums(normalf*normalf))
  d[d<eps] = 1
  normalf = normalf / rbind(t(d),t(d),t(d))
  
  
  # unit normal to the vertex
  for (i in (1:nface)){
    f = F[,i]
    for (j in (1:3)){
      normal[,(f+1)[j]] = normal[,(f+1)[j]]+normalf[,i]
    }
  }
  
  #normalize
  d <- sqrt(colSums(normal*normal))
  d[d<eps] = 1
  normal <- normal / rbind(t(d),t(d),t(d))
  
  # enforce that the normal are outward
  v <- X0 - rbind(colMeans(X0), colMeans(X0) ,colMeans(X0))
  s <- rowSums(v*normal) 
  if (sum(s>0) < sum(s<0)){
    #flip
    normal = -normal
    normalf = -normalf
  }
  
  return (normal)
  
}
