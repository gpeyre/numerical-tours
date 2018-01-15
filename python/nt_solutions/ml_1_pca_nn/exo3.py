# init
I = random.permutation(n)
I = I[0:k]
C = X[I,:]
# Kmeans
niter = 16
it_dist = 0
clf
for it in arange(1,niter+1):
    # NN
    D = pairwise_distances(X,C)
    yb = argmin(D, axis=1)
    # display
    if (it<=3) | (it==niter):
        it_dist = it_dist+1
        subplot(2,2,it_dist)
        for i in arange(0,k):
            I = find(yb==i);
            plot(Z[I,0], Z[I,1], '.')
        CV = (C-mean(X,axis=0)).dot( transpose(V) )
        for i in arange(0,k):
            plot(CV[i,0], CV[i,1], 'ok')
        axis('tight')
        axis('equal')
        axis('off');
    # update centroids
    for l in arange(0,k):
        C[l,:] = mean( X[yb==l,:], axis=0 )
