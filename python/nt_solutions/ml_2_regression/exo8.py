plt.clf
plt.bar(np.arange(1,p+1), abs(xSparse))
plt.bar(np.arange(1,p+1), -abs(xRidge))
plt.legend(('Lasso', 'Ridge'))
