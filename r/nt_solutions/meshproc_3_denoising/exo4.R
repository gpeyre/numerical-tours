ntests <- 15
muList <- seq(from=3, to = 15, by=(15-3)/(ntests-1))/5
errR <- c()

for (i in (1:ntests)){
  mu <- muList[i]
  A <- .sparseDiagonal(x=1, n) + mu*L
  
  for (k in (1:3)){
    Xmu[k,] <- t(pcg(as.matrix(A),X[k,]))
    
  errR <- c(errR,snr(X0,Xmu))
  }
}

plot(errR, type="b", col="blue", xlab = "mu", ylab = "SNR")