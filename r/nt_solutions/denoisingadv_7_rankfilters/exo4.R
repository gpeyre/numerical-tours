

f1 <- f
for (i in 1:4){
  f1 <- opening(f1)
  imageplot(f1, paste("Iteration", i), c(2, 4, i))
}

f1 <- f
for (i in 1:4){
  f1 <- closing(f1)
  imageplot(f1, paste("Iteration", i), c(2, 4, i+4))
}
