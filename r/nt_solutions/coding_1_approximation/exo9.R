




f1 <- fL

for (i in 1:(n%/%w)){
  for (j in 1:(n%/%w)){
    f1[((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w)] <- idct2(fL[((i-1)*w+1):(i*w), ((j-1)*w+1):(j*w)])
  }
}

print(paste("Error |f-f1|/|f| =", norm(as.matrix(f)-f1)/norm(as.matrix(f))))