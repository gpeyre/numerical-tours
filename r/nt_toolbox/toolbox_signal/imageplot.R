




imageplot <- function(f, str='', sbpt=c()){
  ####
  # Use nearest neighbor interpolation for the display.
  ####
  #f <- t( f[1:n,f:n] )
  
  if (length(sbpt) >0){
    if (sbpt[3]==1){
      par(mfrow=sbpt[1:2]) 
    }
  }
  
  if ((length(dim(f)) > 3) && (dim(f)[4] > 1))
  {
    plot(f, main=str, axes=FALSE)
  }
  else {
    f <- as.cimg( t(as.matrix(f)) )
    plot(f, interpolate = FALSE, colorscale = gray, axes = FALSE, main = str)
  }
}
