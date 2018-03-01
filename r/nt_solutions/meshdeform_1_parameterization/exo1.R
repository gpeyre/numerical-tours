options(repr.plot.width=5, repr.plot.height=5)

q=floor(p/4)
t = seq(from=0,to=q-1,by=1)/q
t1 = seq(from=0, to=(p-3*q)-1)/(p-3*q)
Z = rbind(c(t,t*0+1,1-t,t1*0),c(t*0, t, t*0+1,1-t1))
ind = c(seq(from=0,to=p-1),c(0))

plot(x=Z[1,ind+1],y=Z[2,ind+1], type="l",pch=20, col="blue", axes=FALSE, xlab="", ylab="")
points(x=Z[1,ind+1],y=Z[2,ind+1], col="red", pch=20)
