


unique_columns <- function(A){
  sorted_idx <- order(A[1,],A[2,])
  A <- A[,sorted_idx]
  B <- cbind(Matrix(1, nrow=2, ncol=1), t(diff(t(A))))
  R <- pmax(B[1,], B[2,])
  return (A[,which(R!=0)])
}

