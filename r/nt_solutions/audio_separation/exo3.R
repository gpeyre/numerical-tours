options(repr.plot.width=4, repr.plot.height=4)

d = rowSums(P**2)
rho = 0.1
v = d[order(d)]
I = order(d, decreasing=TRUE)

#transformed points
I = I[seq(from=1, to=round(rho*length(I)))]
P1 = P[I,]

#compute Theta
nrow = dim(P1)[1]
Theta = c()
for (i in (1:nrow))
{
  Theta[i] = atan2(P1[i,2], P1[i,1])%%pi
}

nbins = 200
hist = hist(Theta, xlim=c(min(Theta),max(Theta)), col="DarkBlue", main="", breaks=nbins, plot=FALSE)
h = hist$counts/sum(hist$counts)
t = seq(from=pi/200, to = pi, by=(pi-pi/200)/(200-1))
barplot(h, xlab="Theta", col ="DarkBlue", border="DarkBlue"  , tick=TRUE)
