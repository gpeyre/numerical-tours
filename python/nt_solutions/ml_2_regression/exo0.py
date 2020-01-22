niter = 100
Jlist = np.zeros((niter,1))

tau_mult = [.1, .5, 1, 1.8, 1.98];

wopt = np.linalg.solve( X0.transpose().dot(X0), X0.transpose().dot(y0) )
plt.clf

fig, (ax1, ax2) = plt.subplots(2, 1)
for itau in np.arange(0,5):
    tau = tau_mult[itau]/np.linalg.norm(X0,2)**2
    w = np.zeros((p,1))
    for i in np.arange(0,niter):
        Jlist[i] = J(w)
        w = w - tau*GradJ(w)
    #plt.subplot(2,1,1)
    ax1.plot(Jlist)
    ax1.axis('tight')
    plt.title('J(w_k)')
    #plt.subplot(2,1,2)
    e = np.log10( Jlist-J(wopt) +1e-20)
    ax2.plot(e-e[0], label=str(tau_mult[itau]))
    ax2.axis('tight')
    leg = ax2.legend()
    #ax2.legend( str( tau_mult[itau] ) )
    plt.title('$log(J(w_k)-min J)$')
