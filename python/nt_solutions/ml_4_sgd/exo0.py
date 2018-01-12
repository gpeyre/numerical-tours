q=20  # Number of paths computed in parallel.
niter = 1000
x = zeros( (q,niter) )
x[:,0] = random.rand(q,1).flatten()-.5 # Initial conditions.
for i in arange(1,niter):
    u = random.rand(q,1)>.5
    xx = MakeCol(x[:,i-1])
    tau = 1.0/(10+i)
    x1 = xx - tau * ( u*df(xx,1) + (1-u)*df(xx,2) )
    x[:,i] = x1.flatten()
# display
clf
plot(x.transpose(), 'r')
axis('tight')
