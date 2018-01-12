
niter = 20

clf
contourf(tx,ty,F.transpose(),10);

tau_list = array([.5, 1, 1.5, 1.9]) / eta

for itau in arange(0,tau_list.size):
    tau = tau_list[itau]
    x = x0
    X = zeros((2,niter))
    for i in arange(0,niter):
        X[:,i] = x.flatten()
        x = x - tau*Gradf(x)
    plot(X[0,:], X[1,:], '.-');
axis('equal')
