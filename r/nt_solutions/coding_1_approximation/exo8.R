



fL <- matrix(rep(0, length=n*n), c(n, n))

for (i in 1:(n%/%w)){
  for (j in 1:(n%/%w)){
    fL[((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w)] = dct2(as.matrix(f)[((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w)])
  }
}
