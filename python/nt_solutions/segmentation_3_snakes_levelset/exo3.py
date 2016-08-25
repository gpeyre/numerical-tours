plt.figure(figsize=(10,10))
phi = np.copy(phi0) #initialization
eps = np.finfo(float).eps
k = 0

for i in range(1,niter+1):
    g0 = grad(phi, order=2)
    d = np.maximum(eps*np.ones([n,n]), np.sqrt(np.sum(g0**2, 2)))
    g = g0/np.repeat(d[:,:,np.newaxis], 2, 2)
    K = d*div(g[:,:,0], g[:,:,1], order=2)
    phi = phi + tau*K
    if i % int(niter/4.) == 0:
        k = k + 1
        plt.subplot(2, 2, k)
        plot_levelset(phi)