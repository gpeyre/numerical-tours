# init
I = np.random.permutation(n)
I = I[0:k]
C = X[I,:]
# Kmeans
niter = 16
it_dist = 0
plt.clf
for it in np.arange(1,niter+1):
    # NN
    D = pairwise_distances(X,C)
    yb = np.argmin(D, axis=1)
    # display
    if (it<=3) | (it==niter):
        it_dist = it_dist+1
        plt.subplot(2,2,it_dist)
        for i in np.arange(0,k):
            I = np.nonzero(yb==i)[0]
            plt.plot(Z[I,0], Z[I,1], '.')
        CV = (C-np.mean(X,axis=0)).dot( V.transpose() )
        for i in np.arange(0,k):
            plt.plot(CV[i,0], CV[i,1], 'ok')
        plt.axis('tight')
        plt.axis('equal')
        plt.axis('off');
    # update centroids
    for l in np.arange(0,k):
        C[l,:] = np.mean( X[yb==l,:], axis=0 )
