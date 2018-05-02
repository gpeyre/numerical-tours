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
    b = np.ones([N,N,K])
    a = np.copy(b)

    for i in range(niter):

        for k in range(K):
            a[:,:,k] = P[:,:,k]/xi(b[:,:,k])

        q = np.zeros(N)

        for k in range(K):
            q = q + lambd[k] * np.log(np.maximum(1e-19*np.ones(len(b[:,:,k])), b[:,:,k]*xi(a[:,:,k])))

        q = np.exp(q)

        for k in range(K):
            b[:,:,k] = q/xi(a[:,:,k])

    # display
    plt.subplot(m,m,j+1)
    plt.imshow(q)
    plt.axis('off')
