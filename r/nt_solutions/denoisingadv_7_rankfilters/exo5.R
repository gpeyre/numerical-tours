

f1 <- f

for (i in 1:6){
  f1 <- medfilt(f1)
  imageplot(f1, paste("Iteration", i), c(2, 3, i))
}