exo2 <- function(x0, W){
  ####
  # Implement the FM algorithm by iterating these step while the
  # stack |I| is non empty.
  # Display from time to time the front that propagates.
  ####
  n <- dim(W)[1]
  pstart <- t(t(x0))
  DIJSTRA <- perform_dijstra_fm(W, pstart, Inf,'fm', 'sym',n*6)
  D <- DIJSTRA$D ; Dsvg <- DIJSTRA$Dsvg ; Ssvg <- DIJSTRA$Ssvg
  
  par(mfrow=c(2,2))
  for (i in 1:4){
    d <- Dsvg[,,i]
    d[d==Inf] <- 0
    
    d <- as.cimg( t(as.matrix(d)) )
    cmap_jet <- function(v){ return( rgb(v,
                                         (sin(v*2*pi)+1)/2,
                                         (cos(v*2*pi)+1)/2) ) }
    
    plot(d, colourscale=cmap_jet, interpolate = FALSE, axes = FALSE)
    
  }
  
  return(D)
  
}