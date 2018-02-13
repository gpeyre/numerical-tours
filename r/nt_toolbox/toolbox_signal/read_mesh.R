



read_mesh <- function(name){
  ####
  # reading from an OFF file in 3 dimensions, returning X0 (coordinates) and F (faces)
  ####
  ##file <- readLines(name, n = -1L, ok = TRUE, warn= TRUE, encoding = "unknown", skipNul = FALSE)
  
  #check type of file
  if (readLines(name,1) != "OFF") {  
  warning('Wrong type of file, only reads OFF files')
  }
  
  K = lapply(strsplit(readLines(name)," "), as.double)
  
  n_verts = K[[2]][1]
  n_faces = K[[2]][2]
  n_edges = K[[2]][3]
  
  X=matrix(0, nrow=n_verts,ncol=3)
  for (i in (1:n_verts)){
    X[i,]= K[[i+2]]
  }
    
  F=matrix(0, nrow=n_faces, ncol=3)
  for (j in (1:n_faces)){
    F[j,] <- K[[j+n_verts+2]][2:4]
  }
  
  return (list(X0=t(X),F0=t(F)))
}
