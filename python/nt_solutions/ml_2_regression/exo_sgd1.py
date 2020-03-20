xopt = np.linalg.solve( A.transpose().dot(A), A.transpose().dot(y) ) # least square solution
niter = 5000
flist = np.zeros((niter,1))
tau_mult = [.05, .3, .8];
for itau in np.arange(0,len(tau_mult)):
    tauk = tau_mult[itau]/np.linalg.norm(A,2)**2
    x = np.zeros((p,1))
    for i in np.arange(0,niter):
        ik = int( np.floor(np.random.rand()*n) )
        gk = (A[ik,:].dot(x)-y[ik]) * A[ik,:].transpose() # stochastic gradient
        # gk = A.transpose().dot( A.dot(x)-y )  # batch gradient
        x = x - tauk * gk
        flist[i] = f(x)
    plt.plot(flist/f(xopt)-1)
plt.legend(('$\\tau=.05$', '$\\tau=.3$', '$\\tau=.8$'))
