niter = 8000
tau_mult = [.05, .8];
nruns = 10 # number of runs to compute the average performance
for itau in np.arange(0,len(tau_mult)):
    tauk = tau_mult[itau]/np.linalg.norm(A,2)**2
    flist = np.zeros((niter,1))
    for iruns in np.arange(0,nruns):
        x = np.zeros((p,1))
        for i in np.arange(0,niter):
            ik = int( np.floor(np.random.rand()*n) )
            gk = (A[ik,:].dot(x)-y[ik]) * A[ik,:].transpose()
            x = x - tauk * gk
            flist[i] = flist[i] + f(x)
    plt.plot( (flist/nruns)/f(xopt)-1 )
plt.legend(('$\\tau=.05$', '$\\tau=.8$'))
