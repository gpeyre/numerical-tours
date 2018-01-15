



for (i in 1:n){
  D <- pmin(D, array(D[,i],c(n,n)) + t(array(D[i,],c(n,n))))
}