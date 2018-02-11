delta <- 5
sel = c(seq(delta,p-1,1),seq(0,delta-1))

R <- Matrix(0, nrow=2, ncol=n)
R[,(B+1)[sel+1]]=Z

Y <- Matrix(0, nrow=2, ncol=n)
Y[1,]=solve(L1,R[1,])
Y[2,]=solve(L1,R[2,])

plot_mesh(as.matrix(rbind(Y, matrix(0,nrow=1, ncol=n))), F)