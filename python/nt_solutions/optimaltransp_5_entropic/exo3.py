from numpy import linalg

b = np.ones(N)
niter = 2000
Err_p = []
Err_q = []

for i in range(niter):
    a = p/(np.dot(xi,b))
    Err_q = Err_q + [linalg.norm(b*np.dot(xi,a) - q)/linalg.norm(q)]
    b = q /(np.dot(np.transpose(xi),a))
    Err_p = Err_p + [linalg.norm(a*np.dot(np.transpose(xi),b) - p)/linalg.norm(p)]

plt.figure(figsize = (10,7))

plt.subplot(2,1,1)
plt.title("$||\pi -p||$")
plt.plot(np.log(np.asarray(Err_p)), linewidth = 2)

plt.subplot(2,1,2)
plt.title("$||\pi^T -q||$")
plt.plot(np.log(np.asarray(Err_q)), linewidth = 2)

plt.show()                                 