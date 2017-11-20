



tlist <- seq(0.3, 0.9, length=4)

for (i in 1:length(tlist)){
  Th <- tlist[i]
  imageplot(clamp(ThreshBlock(f, Th)),
            paste("T =", round(Th,1)),
            c(2, 2, i))
}