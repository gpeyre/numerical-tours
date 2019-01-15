from numpy import linalg

niter = 1600
v = np.ones([N,N,R])
u = np.copy(v)
Err_q = np.zeros(niter)

for i in range(niter):

    for k in range(R):
        Err_q[i] = Err_q[i] + linalg.norm(u[:,:,k]*K(v[:,:,k]) - A[:,:,k], 1)
        u[:,:,k] = A[:,:,k]/K(v[:,:,k])

    b = np.zeros(N)
    for k in range(R):
        b = b + lambd[k] * np.log(np.maximum(1e-19*np.ones(len(v[:,:,k])), v[:,:,k]*K(u[:,:,k])))
    b = np.exp(b)

    for k in range(R):
        v[:,:,k] = b/K(u[:,:,k])

plt.figure(figsize=(7,5))
plt.plot(np.log(Err_q),linewidth = 2)
plt.show()
