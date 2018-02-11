



cR <- as.vector(abs(fF))[order(-as.vector(abs(fF)))]
plot(log10(cR), type="l", col="blue", lwd=2, xlim=c(0,n**2))
