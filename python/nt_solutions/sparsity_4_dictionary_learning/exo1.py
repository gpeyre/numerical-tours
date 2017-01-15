'''
Sparse code data Y using dictionary D using a forward backward iterative
scheme (projected block coordinate gradient descent).
'''
n_iter = 100
gamma = 1.6/np.linalg.norm(D, ord=2)**2
E = np.zeros((n_iter, m))
X = np.zeros((p, m))
for i in range(n_iter):
    R = np.dot(D, X) - Y
    E[i, :] = np.sum(R**2, axis=0)
    X = projX(X - gamma * np.dot(D.T, R), k)

sel = range(5)
plt.plot(np.log10(E[:len(E)//2, sel] - np.min(E[:, sel], axis=0)))
plt.axis('tight')
plt.title('$log_{10}(J(x_j) - J(x_j^*))$')
plt.show()
