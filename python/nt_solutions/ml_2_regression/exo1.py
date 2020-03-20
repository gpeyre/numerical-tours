q = 50
lmax = np.linalg.norm(A,2)**2
lambda_list = lmax*np.linspace(.3,1e-3,q)
X = np.zeros( (p,q) )
E = np.zeros( (q,1) )
for i in np.arange(0,q):
    Lambda = lambda_list[i]
    x = np.linalg.solve( A.transpose().dot(A) + Lambda*np.eye(p), A.transpose().dot(y) )
    X[:,i] = x.flatten() # bookkeeping
    E[i] = np.linalg.norm(A1.dot(x)-y1) / np.linalg.norm(y1)
# find optimal lambda
i = E.argmin()
lambda0 = lambda_list[i]
xRidge = X[:,i]
print( 'Ridge: ' + str(E.min()*100) + '%' )
# Display error evolution.
plt.clf
plt.plot(lambda_list/lmax, E)
plt.plot( [lambda0/lmax,lambda0/lmax], [E.min(), E.max()], 'r--')
plt.axis('tight')
plt.xlabel('$\lambda/|X|^2$')
plt.ylabel('$E$')
