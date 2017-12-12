

beta_list <- seq(0, 1, length=6)

for (i in 1:length(beta_list)){
  beta_c <- beta_list[i]
  imageplot(phi(f, beta_c), paste("Beta =", beta_c), c(2, 3, i))
}