

mu_list <- seq(0.5,6, length=6)

for (i in (1:length(mu_list))){
  mu <- mu_list[i]
  imageplot( denoise(as.matrix(y),mu), paste("mu =", mu), c(2,3,i) );
}