fig = plt.figure(figsize=(15,15))
niter = 150;
stress = []
Xstress = np.copy(X)
ndisp = [1, 5, 10, min(niter,100), np.float("Inf")]
k = 0

for i in range(niter):

    if ndisp[k] == i:
        ax = fig.add_subplot(2, 2, k+1, projection="3d")
        
        #swiss roll
        ax.scatter(Xstress[0,:], Xstress[1,:], Xstress[2,:], c=plt.cm.jet((X[0,:]**2+X[2,:]**2)/100), s=ms, lw=0, alpha=1)
        
        #graph
        I,J,V = sparse.find(A)
        xx = np.vstack((Xstress[0,I],Xstress[0,J]))
        yy = np.vstack((Xstress[1,I],Xstress[1,J]))
        zz = np.vstack((Xstress[2,I],Xstress[2,J]))
        
        for i in range(len(I)):
            ax.plot(xx[:,i], yy[:,i], zz[:,i], color="black")
        
        #params
        ax.axis("off")
        ax.set_xlim(np.min(Xstress[0,:]),np.max(Xstress[0,:]))
        ax.set_ylim(np.min(Xstress[1,:]),np.max(Xstress[1,:]))
        ax.set_zlim(np.min(Xstress[2,:]),np.max(Xstress[2,:]))
        ax.view_init(elev=el, azim=az)
        
        k += 1

    # Compute the distance matrix.
    D1 = np.repeat(np.sum(Xstress**2, 0)[:,np.newaxis], n, 1)
    D1 = np.sqrt(D1 + np.transpose(D1) - 2*np.dot(np.transpose(Xstress), Xstress))
    
    # Compute the scaling matrix.
    B = -D/np.maximum(D1,1e-10*np.ones(np.shape(D1)))
    B = B - np.diag(np.sum(B,0))
    
    # update
    Xstress = np.transpose(np.dot(B, np.transpose(Xstress)))/n
    # Xstress = Xstress-repmat(mean(Xstress,2), [1 n]);
    # record stress
    stress = stress + [np.sqrt(np.sum(abs(D-D1)**2)/n**2)]

plt.show()