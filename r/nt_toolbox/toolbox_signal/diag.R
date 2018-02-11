



diag <- function(X)
{
  ####
  # Extract the diagonal of a matrix and returns a vector
  ####
  out = rep(0, dim(X)[1])
  for (i in 1:length(out))
  {
    out[i] = X[i,i]
  }
  return(out)
  
}

