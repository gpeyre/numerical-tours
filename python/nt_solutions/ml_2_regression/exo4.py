q = 200
lambda_list = lmax*np.linspace(.6,1e-3,q)
W = np.zeros( (p,q) )
E = np.zeros( (q,1) )
w = np.zeros((p,1))
niter = 500
for iq in np.arange(0,q):
    Lambda = lambda_list[iq]
    # ISTA #
    for i in np.arange(0,niter):
        w = ISTA(w,Lambda,tau)
    W[:,iq] = w.flatten() # bookkeeping
    E[iq] = np.linalg.norm(X1.dot(w)-y1) / np.linalg.norm(y1)
# find optimal Lambda
i = E.argmin()
lambda0 = lambda_list[i]
wSparse = W[:,i]
print( 'Lasso: ' + str(E.min()*100) + '%' )
# Display error evolution.
plt.clf
plt.plot(lambda_list/lmax, E)
plt.plot( [lambda0/lmax,lambda0/lmax], [E.min(),E.max()], 'r--')
plt.axis('tight')
plt.xlabel('$\lambda/|X^* y|_\infty$')
plt.ylabel('$E$')
