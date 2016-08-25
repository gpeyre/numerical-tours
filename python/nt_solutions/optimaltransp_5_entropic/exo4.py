from numpy import linalg

niter = 800
b = np.ones([N,N,K])
a = np.copy(b)
Err_q = np.zeros(niter)

for i in range(niter):
    
    for k in range(K):
        Err_q[i] = Err_q[i] + linalg.norm(a[:,:,k]*xi(b[:,:,k]) - P[:,:,k])/linalg.norm(P[:,:,k])
        a[:,:,k] = P[:,:,k]/xi(b[:,:,k])
        
    q = np.zeros(N)
    
    for k in range(K):
        q = q + lambd[k] * np.log10(np.maximum(1e-19*np.ones(len(b[:,:,k])), b[:,:,k]*xi(a[:,:,k])))

    q = np.exp(q)

    for k in range(K):
        b[:,:,k] = q/xi(a[:,:,k])

plt.figure(figsize=(7,5))
plt.plot(np.log10(Err_q))
plt.show()