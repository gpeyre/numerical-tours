from numpy import linalg

epsilon = .01;
K = np.exp(-C/epsilon)

v = np.ones(N[1])
niter = 1000
Err_p = []
Err_q = []

for i in range(niter):
    # sinkhorn step 1
    u = a / (np.dot(K,v))
    # error computation
    r = v*np.dot(np.transpose(K),u)
    Err_q = Err_q + [linalg.norm(r - b, 1)]
    # sinkhorn step 2
    v = b /(np.dot(np.transpose(K),u))
    s = u*np.dot(K,v)
    Err_p = Err_p + [linalg.norm(s - a,1)]

plt.figure(figsize = (10,7))

plt.subplot(2,1,1)
plt.title("$||P 1 - a||$")
plt.plot(np.log(np.asarray(Err_p) + 1e-5), linewidth = 2)

plt.subplot(2,1,2)
plt.title("$||P^T 1 - b||$")
plt.plot(np.log(np.asarray(Err_q) + 1e-5), linewidth = 2)

plt.show()
