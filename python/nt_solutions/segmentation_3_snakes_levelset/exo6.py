plt.figure(figsize=(10,10))
phi = np.copy(phi0)
k = 0
gW = grad(W, order=2)

for i in range(1,niter+1):
    gD = grad(phi, order=2)
    d = np.maximum(eps*np.ones([n,n]), np.sqrt(np.sum(gD**2, 2)))
    g = gD/np.repeat(d[:,:,np.newaxis], 2, 2)
    G = W*d*div(g[:,:,0], g[:,:,1], order=2) + np.sum(gW*gD,2)
    phi = phi + tau*G
    if i % 30 == 0:
        phi = perform_redistancing(phi)
    if i % int(niter/4.) == 0:
        k = k + 1
        plt.subplot(2, 2, k)
        plot_levelset(phi,0,f0)