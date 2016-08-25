plt.figure(figsize=(20,15))
X1 = X
err = [pnoisy]

for i in range(1,13):
    X1 = tW.dot(np.transpose(X1)).transpose()
    err = err + [snr(X0, X1)]
    
    if i%2 == 0:
        plot_mesh(X1, F, subplot=[2,3,i//2])
        
    if err[-1] > max(err[:-1]):
        Xbest = X1