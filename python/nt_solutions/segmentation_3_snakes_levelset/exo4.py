Y,X = np.meshgrid(np.arange(1,n+1), np.arange(1,n+1))
r = n/3.
c = np.asarray([n,n])/2
phi0 = np.maximum(abs(X-c[0]), abs(Y-c[1])) - r