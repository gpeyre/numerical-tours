Layers = []
Layers.append([p, k]) # vanilla logitic regression (convex)
Layers.append([p, 8, k]) # single hidden layer
Layers.append([p, 3, 4, k]) # 4 hidden layers
Layers.append([p, 4, 5, 8, 4, k]) # 4 hidden layers
tau_list = np.array([.01/10, .01/10, .01/30, .01/40])
plt.clf
for il in np.arange(0,np.size(Layers)):
	D = Layers[il]
	# layers
	R = np.size(D)-1
	A = []
	b = []
	for r in np.arange(0,R):
	    A.append(np.random.randn(D[r+1],D[r]))
	    b.append(np.random.randn(D[r+1],1))
    # descent
	tau = tau_list[il]
	L = [] #
	niter = 12000
	L = np.zeros((niter,1))
	for it in np.arange(0,niter):
	    [L[it],gA,gb] = ForwardBackwardNN(A,b)
	    for r in np.arange(0,R):
	        A[r] = A[r] - tau*gA[r]
	        b[r] = b[r] - tau*gb[r]
    # probability
	V = ForwardNN(A,b,Z,R)
	U = np.reshape(SM(V[R].transpose()), [q,q,k] )
	# same color
	R = np.zeros((q,q,3))
	for i in np.arange(0,k):
		for a in np.arange(0,3):
			R[:,:,a] = R[:,:,a] + U[:,:,i]/np.max( U[:,:,i].flatten() ) * col[a,i]
	# display
	plt.subplot(2,2,il+1)
	M=1
	plt.clf
	plt.imshow(R.transpose((1, 0, 2)), origin="lower", extent=[-M,M,-M,M])
	for i in np.arange(0,k):
	    I = find(y==i+1)
	    plt.plot(x[0,I], x[1,I], '.', color=col[:,i]*.8)
	plt.axis('off');
    # title(['D=' num2str(D)])
