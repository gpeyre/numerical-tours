
exo3 <- function(x0, W){
  ####
  # Compute the distance map to these starting point using the FM algorithm.
  ####
  n <- dim(W)[1]
  pstart <- t(t(x0))
  DIJSTRA <- perform_dijstra_fm(W, pstart, Inf,'fm', 'sym',n*6)
  D <- DIJSTRA$D ; Dsvg <- DIJSTRA$Dsvg ; Ssvg <- DIJSTRA$Ssvg
  # display
  k <- 8
  displ <- function(D){ cos(2*pi*k*D/max(D) ) }
  
  cmap_jet <- function(v){ return( rgb(v,
                                       (sin(v*2*pi)+1)/2,
                                       (cos(v*2*pi)+1)/2) ) }
  
  plot(as.cimg(displ(D)), colourscale=cmap_jet, interpolate = FALSE, axes = FALSE)
  
  return(D)
  
}