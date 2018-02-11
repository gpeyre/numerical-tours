sigma_list = np.array( [.1, .5, 1, 4] )
niter = 4000
plt.clf
for io in np.arange(0, np.size(sigma_list)):
    sigma = sigma_list[io]
    # grad descent
    K = kappa(X,X,sigma)
    tau = .5
    if io==4:
        tau = .05
    h = np.zeros((n,1))
    for i in np.arange(0,niter):
        h = h - tau * nablaF(h,K,y)
    # evaluate on a grid
    K1 = kappa(G,X,sigma)
    Theta = theta( K1.dot(h) )
    Theta = Theta.reshape((q,q))
    # Display the classification probability.
    plt.subplot(2,2,io+1)
    plt.imshow(Theta.transpose(), origin="lower",  extent=[-tmax, tmax, -tmax, tmax])
    plt.plot(X[I,0], X[I,1], '.')
    plt.plot(X[J,0], X[J,1], '.')
    plt.axis('equal')
    plt.axis('off')
    plt.title('$\sigma=' + str(sigma) + '$')
