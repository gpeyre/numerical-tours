



read_mesh <- function(name){
  ####
  # reading from an OFF file in 3 dimensions, returning X0 (coordinates) and F (faces)
  ####
  ##file <- readLines(name, n = -1L, ok = TRUE, warn= TRUE, encoding = "unknown", skipNul = FALSE)
  
  #check type of file
  if (readLines(name,1) != "OFF") {  
  warning('Wrong type of file, only reads OFF files')
  }
  
  n_verts = as.numeric(unlist(strsplit(readLines(name)[2],split=" "))[[1]])
  n_faces = as.numeric(unlist(strsplit(readLines(name)[2],split=" "))[[2]])
  n_edges = as.numeric(unlist(strsplit(readLines(name)[2],split=" "))[[3]])
  
  X=c()
  for (i in (3:(n_verts+2))){
    X <- rbind(X, c(as.numeric(unlist(strsplit(readLines(name)[i],split=" "))[[1]]), 
                       as.numeric(unlist(strsplit(readLines(name)[i],split=" "))[[2]]),
                as.numeric(unlist(strsplit(readLines(name)[i],split=" "))[[3]])))
  }
    
  F=c()
  for (j in (n_verts+3):(n_verts+n_faces+2)){
    F <- rbind(F, c(as.numeric(unlist(strsplit(readLines(name)[j],split=" "))[[2]]), 
                    as.numeric(unlist(strsplit(readLines(name)[j],split=" "))[[3]]),
                    as.numeric(unlist(strsplit(readLines(name)[j],split=" "))[[4]])))
  }
  
  return (list(X0=t(X),F0=t(F)))
}
