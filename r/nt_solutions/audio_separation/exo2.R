sel = sample(n)

sel = sel[1:npts]

plot(P[sel,1], P[sel,2], type="p", pch=19, col = "blue", cex=0.05, main="Time domain", xlab="", ylab="", xlim=c(-5,5), ylim=c(-5,5))