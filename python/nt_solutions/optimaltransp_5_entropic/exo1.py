from numpy import linalg

b = np.ones(N[1])
niter = 300
Err_p = []
Err_q = []

for i in range(niter):
    a = p/(np.dot(xi,b))
    Err_q = Err_q + [linalg.norm(b*np.dot(np.transpose(xi),a) - q)/linalg.norm(q)]
    b = q /(np.dot(np.transpose(xi),a))
    Err_p = Err_p + [linalg.norm(a*np.dot(xi,b) - p)/linalg.norm(p)]

plt.figure(figsize = (10,7))

plt.subplot(2,1,1)
plt.title("$||\pi -p||$")
plt.plot(np.log10(np.asarray(Err_p) + 1e-5), linewidth = 2)

plt.subplot(2,1,2)
plt.title("$||\pi^T -q||$")
plt.plot(np.log10(np.asarray(Err_q) + 1e-5), linewidth = 2)

plt.show()  