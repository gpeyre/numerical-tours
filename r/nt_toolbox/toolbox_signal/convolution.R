convolution = function(X, A)
{
  X.pad = matrix(0, ncol=(NCOL(X) + NCOL(A) - 1), nrow=(NROW(X) + NROW(A) - 1))
  X.pad[1:NROW(X), 1:NCOL(X)] = X
  A.pad = matrix(0, ncol=(NCOL(X) + NCOL(A) - 1), nrow=(NROW(X) + NROW(A) - 1))
  A.pad[1:NROW(A), 1:NCOL(A)] = A

  X.fft = fft(X.pad)
  A.fft = fft(A.pad)
  M = fft(X.fft * A.fft, inverse = TRUE) / length(X.fft)

  N.row = (floor(NROW(A)/2) + (1:NROW(X)))
  N.col = (floor(NCOL(A)/2) + (1:NCOL(X)))

  XC = M[N.row, N.col]
  return(Re(XC))
}