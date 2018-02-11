compute_boundary <- function(face){
  ###   compute_boundary - compute the vertices on the boundary of a 3D mesh
  
  ###   boundary=compute_boundary(face);
  

  #if dim(face)[1] < dim(face)[2]
    #face<-t(face)
  
  # compute edges (i,j) that are adjacent to only 1 face
  A=compute_edge_face_ring(face)
  i<-which(A!=0, arr.ind=TRUE)[,2]
  j<-which(A!=0, arr.ind=TRUE)[,1]
  v<-t(A)[A!=0]
    
  
  i = i[v==-1]
  j= j[v==-1]
  
  
  # build the boundary by traversing the edges
  boundary <- c(i[1])
  i<-i[-1]
  j<-j[-1]
  
  while (length(i)>0){
    b<-boundary[length(boundary)]
    I<-which(i==0)
      
      if (length(I) == 0){
          I <- which(j==b)
          
          if (length(I)==0){
            warning('Problem with boundary')
          }
        
          boundary <- c(boundary,i[I][1])
      }
      else{
          boundary <- c(boundary,j[I][1])
      }
          i<-i[-I]
          j<-j[-I]
          
    }
  return (boundary-1) 
}