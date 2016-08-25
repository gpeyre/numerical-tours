plt.figure(figsize=(15,15))
Xt = np.copy(X)
k = 0
sob = []
err = []

for i in range(1,niter+1):
    # step
    Xt = Xt - tau*(tL.dot(np.transpose(Xt))).transpose()
    # error
    err = err + [snr(X0, Xt)]
    if i % np.floor(niter/4) == 0:
        k = k + 1
        plot_mesh(Xt, F, subplot=[2,2,k], title="T = %i" %(Tmax*k/4))