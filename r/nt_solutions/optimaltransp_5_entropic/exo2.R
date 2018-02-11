



glist <- c(.1,.01,.005,.001)
niter <- 300

for (k in 1:length(glist)){
  gamma <- glist[k]
  xi <- exp(-C/gamma)
  b <- rep(1, N[2])
  
  for (i in 1:niter){
    a <- p / (xi %*% b)
    b <- q / (t(xi) %*% a)
  }
  
  a_diag <- array(0, c(length(a), length(a)))
  diag(a_diag) <- a
  
  b_diag <- array(0, c(length(b), length(b)))
  diag(b_diag) <- b
  
  Pi <- (a_diag %*% xi) %*% b_diag
  
  imageplot(clamp(Pi,0,min(1/N)*.3), paste("gamma=", gamma), c(2,2,k))
  
}