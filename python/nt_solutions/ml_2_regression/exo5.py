plt.clf
for i in np.arange(0,p):
    plt.plot(lambda_list/lmax, W[i,:], label=class_names[0][i])
plt.plot( [lambda0/lmax,lambda0/lmax], [W.flatten().min(),W.flatten().max()], 'r--')
plt.axis('tight')
plt.xlabel('$\lambda/|X^* y|_\infty$')
plt.ylabel('$w_i$')
plt.legend()
