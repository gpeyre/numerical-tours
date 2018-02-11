def ForwardBackwardNN(A,b):
    ## forward pass
    X = ForwardNN(A,b,x,R)
    L = lossF(X[R],Y)
    [gA,gb] = BackwardNN(A,b,X,R)
    return [L,gA,gb]

tau = .01/5
tau = .01/10
if R==5:
    tau = .01/50
if R==3:
    tau = .01/80
niter = 8000
L = np.zeros((niter,1))
for it in np.arange(0,niter):
    [L[it],gA,gb] = ForwardBackwardNN(A,b)
    for r in np.arange(0,R):
        A[r] = A[r] - tau*gA[r]
        b[r] = b[r] - tau*gb[r]
plt.clf
plt.plot(L)
plt.xlabel('iter')
plt.ylabel('$L$')
plt.axis('tight')
