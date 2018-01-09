



slist <- c(1, 2, 4, 6)

for (i in 1:length(slist)){
  sigma <- slist[i]
  d <- sqrt(apply(nabla(blur(as.matrix(f0), sigma))**2, c(1,2), sum))
  t <- max(d)*1./5
  imageplot(d > t, paste("sigma =", format(sigma, digits=1)), c(2,2,i))
}