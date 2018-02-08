



slist <- c(4, 6, 10, 15)

for (i in 1:length(slist)){
  sigma = slist[i]
  plot_levelset(delta(blur(f0, sigma)) , f0, paste("sigma =", sigma), sbpt=c(2,2,i))
}