sigma_list = [.05 .1 .5 1 5 10]

for i=1:length(sigma_list)
    sigma = sigma_list[i]
    kappa = (X,Z) -> exp( -distmat(X,Z)/(2*sigma^2) )
    # Regressor.
    h = (kappa(X,X)+lambda*eye(n))\y
    Y = x -> kappa(x,X)*h
    # Eval on the grid
    yn = reshape(Y(Xn),q,q);
    #
    subplot(2,3,i)
    imshow(yn, extent=[0,1, 0,1], cmap = get_cmap("jet"))
    axis("image"); axis("off"); 
    title(string("\sigma=",sigma));
end
