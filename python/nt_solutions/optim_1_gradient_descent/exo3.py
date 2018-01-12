
niter = 20

x = x0

X = zeros((2,niter))
for i in arange(0,niter):
    X[:,i] = x.flatten()
    # optimal step size
    tau = (x[0]**2+eta**2*x[1]**2)/(x[0]**2+eta**3*x[1]**2)
    x = x - tau*Gradf(x)


clf
contourf(tx,ty,F.transpose(),10);
plot(X[0,:], X[1,:], 'k.-');
axis('equal')
