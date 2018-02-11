# Defining the L1 norm of a vector


norm <- function(v){
  ####
  # Euclidian norm of a vector
  ####
  return(sqrt(sum(v**2)))
}



l1_norm = function(x){sum(abs(x))}


