a = 5
b = 3.5
tx = np.linspace(-a, a, num=q)
ty = np.linspace(-b, b, num=q)
[B, A] = np.meshgrid(ty, tx)
G = np.vstack([A.flatten(), B.flatten()]).transpose()
##
offs = np.array([0.3, 1, 3, 5]);
niter = 10000
plt.clf
for io in np.arange(0, np.size(offs)):
    # generate data
    omega = offs[io]*np.array([1, .5])
    X = np.vstack((np.random.randn(n1,2)-np.ones([n1,1])*omega/2, np.random.randn(n1,2)+np.ones([n1,1])*omega/2))
    # run gradient descent
    w = np.zeros((p+1, 1))  # initialization
    for i in np.arange(0, niter):
        w = w - tau * nablaE(w, AddBias(X), y)
    # display
    Theta = theta(AddBias(G).dot(w))
    Theta = Theta.reshape((q,q))
    plt.subplot(2, 2, io+1)
    plt.imshow(Theta.transpose(), origin="lower",  extent=[-a, a, -b, b])
    plt.axis('equal')
    plt.plot(X[I,0], X[I,1], '.')
    plt.plot(X[J,0], X[J,1], '.')
    plt.axis('off')
    plt.axis([-a, a, -b, b])
