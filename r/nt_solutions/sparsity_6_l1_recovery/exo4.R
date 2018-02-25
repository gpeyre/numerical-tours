klist = c(10, 30, 50)
P = 200
ntrials = 200
tmin = 0
tmax = 2.5
q = 50
t = seq(tmin, tmax, length=q)
t1 = seq(tmin, tmax, length=1000)
dt = (tmax-tmin)/q

options(repr.plot.width=6, repr.plot.height=3)

for (j in 1:length(klist))
{
    k = klist[j]
    
    # simulation    
    v = c()
    for (i in 1:ntrials)
    {
        v = c(v, svd(randn(P, k)/sqrt(P))$d **2)
    }
    #plt.figure(figsize = (10,10))
    #plt.subplot(len(klist),1,j+1)
    h = hist(v, breaks=t, plot=FALSE)$counts
    h = h/sum(h)/dt
    barplot(h, main=paste("P=", P, ", k=", k), col=4, xlab='', ylab='')
    
    # theoritical law
    beta = k/P
    a = (1 - sqrt(beta))**2
    b = (1 + sqrt(beta))**2
    z = sqrt(pmax(t1-a, rep(0, length(t1))) * pmax(b - t1, rep(0, length(t1))))/(2*pi*beta*t1)
    par(new=TRUE)
    plot(t1, z, yaxt='n', type="l", col=2, xlab='', ylab='')
    
}
