a = 5
b = 3.5
tx = linspace(-a, a, num=q)
ty = linspace(-b, b, num=q)
[B, A] = meshgrid(ty, tx)
G = transpose(np.vstack([A.flatten(), B.flatten()]))
##
offs = array([0.3, 1, 3, 5]);
niter = 10000
clf
for io in range(0, size(offs)):
    # generate data
    omega = offs[io]*array([1, .5])
    X = np.vstack((np.random.randn(n1,2)-np.ones([n1,1])*omega/2, np.random.randn(n1,2)+np.ones([n1,1])*omega/2))
    # run gradient descent
    w = zeros((p+1, 1))  # initialization
    for i in range(0, niter):
        w = w - tau * nablaE(w, AddBias(X), y)
    # display
    Theta = theta(AddBias(G).dot(w))
    Theta = Theta.reshape((q,q))
    subplot(2, 2, io+1)
    imshow(transpose(Theta), origin="lower",  extent=[-a, a, -b, b])
    axis('equal')
    plot(X[I,0], X[I,1], '.')
    plot(X[J,0], X[J,1], '.')
    axis('off')
    axis([-a, a, -b, b])
