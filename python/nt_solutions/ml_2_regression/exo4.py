q = 200
lambda_list = lmax*linspace(.6,1e-3,q)
W = zeros( (p,q) )
E = zeros( (q,1) )
w = zeros((p,1))
niter = 500
for iq in arange(0,q):
    Lambda = lambda_list[iq]
    # ISTA #
    for i in arange(0,niter):
        w = ISTA(w,Lambda,tau)
    W[:,iq] = w.flatten() # bookkeeping
    E[iq] = norm(X1.dot(w)-y1) / norm(y1)
# find optimal Lambda
i = E.argmin()
lambda0 = lambda_list[i]
# Display error evolution.
clf
plot(lambda_list/lmax, E)
plot( [lambda0/lmax,lambda0/lmax], [E.min(),E.max()], 'r--')
axis('tight')
xlabel('$\lambda/|X^* y|_\infty$')
ylabel('$E$')
