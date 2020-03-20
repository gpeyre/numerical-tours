plt.clf
for i in np.arange(0,p):
    plt.plot(lambda_list/lmax, X[i,:], label=class_names[0][i])
plt.plot( [lambda0/lmax,lambda0/lmax], [X.flatten().min(),X.flatten().max()], 'r--')
plt.axis('tight')
plt.xlabel('$\lambda/|A^* y|_\infty$')
plt.ylabel('$x_i$')
plt.legend()
