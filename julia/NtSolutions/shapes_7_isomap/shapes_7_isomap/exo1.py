for i in range(n):
    D = np.minimum(D, np.tile(D[:,i][:,np.newaxis], (1,n)) + np.tile(D[i,:], (n,1)))