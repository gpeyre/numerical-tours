q=20  # Number of paths computed in parallel.
niter = 1000
x = np.zeros( (q,niter) )
x[:,0] = np.random.rand(q,1).flatten()-.5 # Initial conditions.
for i in np.arange(1,niter):
    u = np.random.rand(q,1)>.5
    xx = MakeCol(x[:,i-1])
    tau = 1.0/(10+i)
    x1 = xx - tau * ( u*df(xx,1) + (1-u)*df(xx,2) )
    x[:,i] = x1.flatten()
# display
plt.clf
plt.plot(x.transpose(), 'r')
plt.axis('tight')
