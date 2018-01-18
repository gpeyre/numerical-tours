plt.clf
plt.bar(np.arange(1,p+1), abs(wSparse))
plt.bar(np.arange(1,p+1), -abs(wRidge))
plt.legend(('Lasso', 'Ridge'))
