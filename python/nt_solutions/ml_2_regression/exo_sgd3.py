niter = 20000
flist = np.zeros((niter,1))
x = np.zeros((p,1))
for i in np.arange(0,niter):
    tauk = 1/np.linalg.norm(A,2)**2 * 1/(1+i/10)
    ik = int( np.floor(np.random.rand()*n) )
    gk = (A[ik,:].dot(x)-y[ik]) * A[ik,:].transpose()
    x = x - tauk * gk
    flist[i] = f(x)
plt.plot(flist/f(xopt)-1 );
print( flist[-1]/f(xopt)-1 )
