sigma_list = [.1 .5 1 4]
niter = 4000

for is=1:length(sigma_list)
    sigma = sigma_list[is]
    # grad descent
    K = kappa(X,X,sigma)
    Flist = []
    tau = .5
    if is==4
        tau = .05
    end
    h = zeros(n,1)
    for i=1:niter
        h = h - tau * nablaF(h,K,y)
        append!(Flist, F(h,K,y))
    end
    # evaluate on a grid
    Theta = reshape( theta(kappa(G,X,sigma)*h) , q,q)
    # Display the classification probability.
    subplot(2,2,is)
    imshow(Theta'[:,end:-1:1], extent=[-tmax, tmax, -tmax, tmax], cmap = get_cmap("jet")); #cmap parula
    plot_multiclasses(X,y,ms=5); axis("off")
    clim(0, 1)
    title(string("\sigma=", sigma))
end
