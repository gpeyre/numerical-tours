sigma_list = array( [.1, .5, 1, 4] )
niter = 4000
clf
for io in arange(0, size(sigma_list)):
    sigma = sigma_list[io]
    # grad descent
    K = kappa(X,X,sigma)
    tau = .5
    if io==4:
        tau = .05
    h = zeros((n,1))
    for i in arange(0,niter):
        h = h - tau * nablaF(h,K,y)
    # evaluate on a grid
    K1 = kappa(G,X,sigma)
    Theta = theta( K1.dot(h) )
    Theta = Theta.reshape((q,q))
    # Display the classification probability.
    subplot(2,2,io+1)
    imshow(transpose(Theta), origin="lower",  extent=[-tmax, tmax, -tmax, tmax])
    plot(X[I,0], X[I,1], '.')
    plot(X[J,0], X[J,1], '.')
    axis('equal')
    axis('off')
    title('$\sigma=' + str(sigma) + '$')
