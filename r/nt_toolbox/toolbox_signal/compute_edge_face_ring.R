
compute_edge_face_ring <- function(faces){
  ###   compute_edge_face_ring - compute faces adjacent to each edge
  
  ###   e2f = compute_edge_face_ring(faces);
  
  ###   e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to
  ###   edge (i,j).
  
  n=max(faces)+1
  m=ncol(faces)
  
  X=as.matrix(faces)
  i <- c(X[1,],X[2,],X[3,])
  j <- c(X[2,],X[3,],X[1,])
  s <- c(c(0:(m-1)),c(0:(m-1)),c(0:(m-1)))
  
  # First without duplicate
  tmpt <- (i+(max(i)+1)*j)[!duplicated(i+(max(i)+1)*j)]
  tmp<-sort(tmpt)
  I <- match(tmp,tmpt)-1
  
  # remaining items
  
  J <- setdiff(c(0:(length(s)-1)),I)
  
  # flip the duplicates
  i1 <- c(i[I+1],j[J])
  j1 <- c(j[I+1],i[J])
  s <- c(s[I+1],s[J])
  
  
  # remove doublons
  tmpt <- (i1+(max(i1)+1)*j1)[!duplicated(i1+(max(i1)+1)*j1)]
  tmp<-sort(tmpt)
  I <- match(tmp,tmpt)-1

  i1 <- i1[I+1]
  j1 <- j1[I+1]
  s  <- s[I+1]
  A <- sparseMatrix(i=i1+1, j=j1+1,x=s+1)
  
  # add missing points
  B<- as.matrix(A)
  I<- which(as.vector(B)!=0)
  J<- which(as.vector(t(B))[I]==0)
  C <- replace(as.vector(t(B)),I[J],-1)
  
  
  return (t(matrix(C,nrow=n,ncol=n)))
  
}
