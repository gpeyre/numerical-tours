Layers = []
Layers.append([p, k]) # vanilla logitic regression (convex)
Layers.append([p, 8, k]) # single hidden layer
Layers.append([p, 3, 4, k]) # 4 hidden layers
Layers.append([p, 4, 5, 8, 4, k]) # 4 hidden layers
tau_list = array([.01/10, .01/10, .01/30, .01/40])
clf
for il in arange(0,size(Layers)):
	D = Layers[il]
	# layers
	R = size(D)-1
	A = []
	b = []
	for r in arange(0,R):
	    A.append(randn(D[r+1],D[r]))
	    b.append(randn(D[r+1],1))
    # descent
	tau = tau_list[il]
	L = [] #
	niter = 12000
	L = zeros((niter,1))
	for it in arange(0,niter):
	    [L[it],gA,gb] = ForwardBackwardNN(A,b)
	    for r in arange(0,R):
	        A[r] = A[r] - tau*gA[r]
	        b[r] = b[r] - tau*gb[r]
    # clf plot(1:niter, L, 'LineWidth', 2) axis tight
    # probability
	V = ForwardNN(A,b,Z)
	U = reshape(SM(V[R].transpose()), [q,q,k] )
	# same color
	R = zeros((q,q,3))
	for i in arange(0,k):
		for a in arange(0,3):
			R[:,:,a] = R[:,:,a] + U[:,:,i]/max( U[:,:,i].flatten() ) * col[a,i]
	# display
	subplot(2,2,il+1)
	M=1
	clf
	imshow(R.transpose((1, 0, 2)), origin="lower", extent=[-M,M,-M,M])
	for i in arange(0,k):
	    I = find(y==i+1)
	    plot(x[0,I], x[1,I], '.', color=col[:,i]*.8)
	axis('off');
    # title(['D=' num2str(D)])
