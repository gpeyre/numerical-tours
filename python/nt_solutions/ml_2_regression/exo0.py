niter = 100
flist = np.zeros((niter,1))

tau_mult = [.1, .5, 1, 1.8, 1.98];

xopt = np.linalg.solve( A.transpose().dot(A), A.transpose().dot(y) )
plt.clf

fig, (ax1, ax2) = plt.subplots(2, 1)
for itau in np.arange(0,5):
    tau = tau_mult[itau]/np.linalg.norm(A,2)**2
    x = np.zeros((p,1))
    for i in np.arange(0,niter):
        flist[i] = f(x)
        x = x - tau*Gradf(x)
    #plt.subplot(2,1,1)
    ax1.plot(flist)
    ax1.axis('tight')
    plt.title('f(x_k)')
    #plt.subplot(2,1,2)
    e = np.log10( flist-f(xopt) +1e-20)
    ax2.plot(e-e[0], label=str(tau_mult[itau]))
    ax2.axis('tight')
    leg = ax2.legend()
    #ax2.legend( str( tau_mult[itau] ) )
    plt.title('$log(f(x_k)-min J)$')
