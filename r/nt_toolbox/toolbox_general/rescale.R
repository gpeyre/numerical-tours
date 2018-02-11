




rescale <- function(f, a=0, b=1){
  ####
  # Rescale linearly the dynamic of a vector to fit within a range [a,b]
  ####
  v <- max(f) - min(f)
  g <- (f - min(f))
  if (v>0){
    g <- g/v
  }
  return(a + g*(b-a))
}
