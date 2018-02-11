



names <- c("regular3", "phantom", "lena", "mandrill")
Err_list <- array(rep(0, length=length(11:round(n*n/10))*4), c(4, length(11:round(n*n/10))))


for (i in 1:dim(fList)[3]){
  
  fW <- perform_wavortho_transf(fList[ ,  , i], Jmin, + 1, h)
  cR <- as.vector(abs(fW))[order(-as.vector(abs(fW)))]
  err <- norm(fList[,,i])**2 - cumsum(cR**2)*(norm(fList[,,i])**2/sum(cR**2))
  Err <- err[11:round(n*n/10)]
  
  
  Err_list[i,] <- Err

  }


plot(log10(11:round(n*n/10)), log10(Err_list[1,]/Err_list[1,1]), xlim=c(1,log10(round(n*n/10))), ylim=c(-7,0), type="l", col="blue", lwd=2)
lines(log10(11:round(n*n/10)), log10(Err_list[2,]/Err_list[2,1]), type="l", col="orange", lwd=2)
lines(log10(11:round(n*n/10)), log10(Err_list[3,]/Err_list[3,1]), type="l", col="red", lwd=2)
lines(log10(11:round(n*n/10)), log10(Err_list[4,]/Err_list[4,1]), type="l", col="purple", lwd=2)

title(main="log_10(epsilon^2[M]/ ||f||^2)")
legend(x="topright",
       legend=names,
       col=c("blue", "orange", "red", "purple"),
       inset = 0.05,
       lwd=2,
       cex=1)