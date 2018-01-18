niter = 20

plt.clf
plt.contourf(tx,ty,F.transpose(),10);

tau_list = np.array([.5, 1, 1.5, 1.9]) / eta

for itau in np.arange(0,tau_list.size):
    tau = tau_list[itau]
    x = x0
    X = np.zeros((2,niter))
    for i in np.arange(0,niter):
        X[:,i] = x.flatten()
        x = x - tau*Gradf(x)
    plt.plot(X[0,:], X[1,:], '.-');
plt.axis('equal')
