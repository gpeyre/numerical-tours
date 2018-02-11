s1 = np.hstack((np.arange(2,nbins+1),np.array([nbins-1]))) - 1
s2 = np.hstack((np.array([2]),(np.arange(1,nbins)))) - 1
I = np.where((h[s1] < h ) & (h[s2] < h))
v = np.sort(h[I])
u = np.argsort(h[I])[::-1]
theta1 = t[I[0][u[:3]]]
M1 = np.vstack((np.cos(theta1), np.sin(theta1)))
print("--- M ---")
print(M)
print("--- M1 ---")
print(M1)