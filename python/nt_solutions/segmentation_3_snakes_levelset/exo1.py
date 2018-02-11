r = n/3.
c = n - 10 - np.array([r,r])
phi2 = np.maximum(abs(X-c[0]), abs(Y-c[1])) - r