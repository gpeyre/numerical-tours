niter = 20000
flist = np.zeros((niter,1))
x  = np.zeros((p,1))
x1 = np.zeros((p,1))
for i in np.arange(0,niter):
    tauk = 1/np.linalg.norm(A,2)**2 * 1/(1+np.sqrt(i/10.0))
    ik = int( np.floor(np.random.rand()*n) )
    gk = (A[ik,:].dot(x)-y[ik]) * A[ik,:].transpose()
    x = x - tauk * gk
    x1 = 1/(i+1) * x + i/(i+1) * x1
    flist[i] = f(x1)
plt.plot(flist/f(xopt)-1 );
