q = 50
lmax = norm(X0,2)**2
lambda_list = lmax*linspace(.3,1e-3,q)
W = zeros( (p,q) )
E = zeros( (q,1) )
for i in arange(0,q):
    Lambda = lambda_list[i]
    w = inv(X0.transpose().dot(X0) + Lambda*eye(p)).dot( X0.transpose().dot(y0) )
    W[:,i] = w.flatten() # bookkeeping
    E[i] = norm(X1.dot(w)-y1) / norm(y1)
# find optimal lambda
i = E.argmin()
lambda0 = lambda_list[i]
wRidge = W[:,i]
print( 'Ridge: ' + str(E.min()*100) + '%' )
# Display error evolution.
clf
plot(lambda_list/lmax, E)
plot( [lambda0/lmax,lambda0/lmax], [E.min(), E.max()], 'r--')
axis('tight')
xlabel('$\lambda/|X|^2$')
ylabel('$E$')
