




slist <- c(1,2,5,10)

for (i in 1:length(slist)){
  sigma <- slist[i]
  imageplot(blur(as.matrix(f0), sigma), paste("sigma =", sigma), c(2,2,i))
}