from numpy import linalg

J = np.identity(n) - np.ones([n,n])/n
K = -1/2.*np.dot(J,np.dot(D**2,J))

[val, Xstrain] = sparse.linalg.eigs(K, k=2, which='LR')

Xstrain = Xstrain*np.tile(np.sqrt(val), (n,1))
Xstrain = np.real(np.transpose(Xstrain))

#plot size
plt.figure(figsize = (15,6))

#plot points
plt.scatter(Xstrain[0,:], Xstrain[1,:], ms, c=plt.cm.jet((X[0,:]**2+X[2,:]**2)/100), lw=0, alpha=1)

#plot vertices
I,J,V = sparse.find(A)
xx = np.vstack((Xstrain[0,I], Xstrain[0,J]))
yy = np.vstack((Xstrain[1,I], Xstrain[1,J]))

for i in range(len(I)):
    plt.plot(xx[:,i], yy[:,i], color="black")
    
#params
plt.axis("off")
plt.xlim(np.min(Xstrain[0,:]-1),np.max(Xstrain[0,:])+1)
plt.ylim(np.min(Xstrain[1,:]-1),np.max(Xstrain[1,:])+1)

plt.show()