




rho_list <- c(2,4,8,16)

for (i in 1:length(rho_list)){
  rho <- rho_list[i]
  imageplot(W(f, rho*U), paste("rho =", rho), c(2,2,i))
}