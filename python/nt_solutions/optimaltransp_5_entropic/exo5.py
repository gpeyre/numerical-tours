# Insert your code here.

m = 5
[T,S] = np.meshgrid(np.linspace(0,1,m), np.linspace(0,1,m))
T = np.ravel(T,order="F")
S = np.ravel(S,order="F")
niter = 1000

plt.figure(figsize=(10,5))
for j in range(m**2):
    # weights
    lambd = np.hstack((S[j]*T[j], (1-S[j])*T[j], S[j]*(1-T[j]), (1-S[j])*(1-T[j])))
    # computation
    v = np.ones([N,N,R])
    u = np.copy(v)

    for i in range(niter):
        for k in range(R):
            u[:,:,k] = A[:,:,k]/K(v[:,:,k])
        b = np.zeros(N)
        for k in range(R):
            b = b + lambd[k] * np.log(np.maximum(1e-19*np.ones(len(v[:,:,k])), v[:,:,k]*K(u[:,:,k])))
        b = np.exp(b)
        for k in range(R):
            v[:,:,k] = b/K(u[:,:,k])

    # display
    plt.subplot(m,m,j+1)
    plt.imshow(b)
    plt.axis('off')
