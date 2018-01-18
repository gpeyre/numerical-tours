niter = 20

x = x0

X = np.zeros((2,niter))
for i in np.arange(0,niter):
    X[:,i] = x.flatten()
    # optimal step size
    tau = (x[0]**2+eta**2*x[1]**2)/(x[0]**2+eta**3*x[1]**2)
    x = x - tau*Gradf(x)


plt.clf
plt.contourf(tx,ty,F.transpose(),10);
plt.plot(X[0,:], X[1,:], 'k.-');
plt.axis('equal')
