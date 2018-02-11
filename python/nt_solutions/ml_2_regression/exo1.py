q = 50
lmax = np.linalg.norm(X0,2)**2
lambda_list = lmax*np.linspace(.3,1e-3,q)
W = np.zeros( (p,q) )
E = np.zeros( (q,1) )
for i in np.arange(0,q):
    Lambda = lambda_list[i]
    w = np.linalg.solve( X0.transpose().dot(X0) + Lambda*np.eye(p), X0.transpose().dot(y0) )
    W[:,i] = w.flatten() # bookkeeping
    E[i] = np.linalg.norm(X1.dot(w)-y1) / np.linalg.norm(y1)
# find optimal lambda
i = E.argmin()
lambda0 = lambda_list[i]
wRidge = W[:,i]
print( 'Ridge: ' + str(E.min()*100) + '%' )
# Display error evolution.
plt.clf
plt.plot(lambda_list/lmax, E)
plt.plot( [lambda0/lmax,lambda0/lmax], [E.min(), E.max()], 'r--')
plt.axis('tight')
plt.xlabel('$\lambda/|X|^2$')
plt.ylabel('$E$')
