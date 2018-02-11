



t_list <- max(d)*c(1./4,1./5,1./10,1./20)

for (i in 1:length(t_list)){
  t <- t_list[i]
  imageplot(d > t, paste("t =", format(t, digits=1)) , c(2,2,i))
}