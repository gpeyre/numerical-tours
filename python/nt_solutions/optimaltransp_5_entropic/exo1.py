from numpy import linalg

gamma = .01;
xi = np.exp(-C/gamma)

b = np.ones(N[1])
niter = 1000
Err_p = []
Err_q = []

for i in range(niter):
    # sinkhorn step 1
    a = p/(np.dot(xi,b))
    # error computation
    r = b*np.dot(np.transpose(xi),a)
    Err_q = Err_q + [linalg.norm(r - q, 1)]
    # sinkhorn step 2
    b = q /(np.dot(np.transpose(xi),a))
    s = a*np.dot(xi,b)
    Err_p = Err_p + [linalg.norm(s - p,1)]

plt.figure(figsize = (10,7))

plt.subplot(2,1,1)
plt.title("$||\pi 1 -p||$")
plt.plot(np.log(np.asarray(Err_p) + 1e-5), linewidth = 2)

plt.subplot(2,1,2)
plt.title("$||\pi^T 1 -q||$")
plt.plot(np.log(np.asarray(Err_q) + 1e-5), linewidth = 2)

plt.show()
