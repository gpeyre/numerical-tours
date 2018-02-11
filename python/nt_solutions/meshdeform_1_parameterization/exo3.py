delta = 5
sel = np.hstack((np.arange(delta,p),np.arange(0,delta)))

R = np.zeros([2,n])
R[:,np.asarray(B)[sel]] = Z

Y = np.zeros([2,n])
Y[0,:] = linalg.spsolve(L1, R[0,:])
Y[1,:] = linalg.spsolve(L1, R[1,:])

plt.figure(figsize=(10,10))
plot_mesh(np.vstack((Y,np.zeros(n))),F, lwdt=1, c="lightgrey")