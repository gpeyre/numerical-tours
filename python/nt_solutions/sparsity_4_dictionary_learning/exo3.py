n_iter_learning = 10
n_iter_dico = 50
n_iter_coef = 100
E0 = np.zeros(2*n_iter_learning)
X = np.zeros((p, m))
D = D0.copy()
for i in range(n_iter_learning):
    # --- coefficient update ----
    E = np.zeros((n_iter_coef, m))
    gamma = 1.6/np.linalg.norm(D, ord=2)**2
    for j in range(n_iter_coef):
        R = np.dot(D, X) - Y
        E[j, :] = np.sum(R**2, axis=0)
        X = projX(X - gamma * np.dot(D.T, R), k)
    E0[2*i] = np.linalg.norm(Y - np.dot(D, X))**2
    # --- dictionary update ----
    E = np.zeros(n_iter_dico)
    tau = 1/np.linalg.norm(np.dot(X, X.T), ord=2)
    for j in range(n_iter_dico):
        R = np.dot(D, X) - Y
        E[i] = np.sum(R**2)
        D = projC(D - tau * np.dot(R, X.T))
    E0[2*i+1] = np.linalg.norm(Y - np.dot(D, X))**2

index = list(range(E0.shape[0]))
index_coef = list(range(0, E0.shape[0], 2))
index_dict = list(range(1, E0.shape[0], 2))
plt.plot(np.divide(index, 2), E0, label='$|Y-DX|^2$')
plt.plot(np.divide(index_coef, 2), E0[index_coef],
         '*', markersize=7, label='After coefficient update')
plt.plot(np.divide(index_dict, 2), E0[index_dict],
         'o', markersize=5, label='After dictionary update')
plt.legend(numpoints=1)
plt.axis('tight')
plt.show()
