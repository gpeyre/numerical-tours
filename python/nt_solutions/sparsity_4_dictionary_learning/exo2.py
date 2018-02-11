n_iter = 100
E = np.zeros(n_iter)
tau = 1/np.linalg.norm(np.dot(X, X.T), ord=2)
for i in range(n_iter):
    R = np.dot(D, X) - Y
    E[i] = np.sum(R**2)
    D = projC(D - tau * np.dot(R, X.T))

plt.plot(np.log10(E[:len(E)//2] - np.min(E)))
plt.axis('tight')
plt.show()
