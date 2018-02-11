library(tuneR)

load_sound <- function(file, n0){
  
  x_raw <- readWave(file)@left
  
  n <- length(x_raw)
  
  if (file == "nt_toolbox/data/bird.wav"){
   x_raw <- c(x_raw[6002:12499], x_raw[15001:22499], x_raw[24001:34000]) 
  }
   
  if((n0 < n) & (n0 !=0)){
    x_raw <- x_raw[1:n0]
  }
  return (x_raw/max(x_raw))
}
