sigma_list = [.1 .5 1 4];
niter = 4000;
clf;
for is=1:length(sigma_list)
    sigma = sigma_list(is);
    % grad descent
    K = kappa(X,X,sigma);
    Flist = [];
    tau = .5;
    if is==4
        tau = .05;
    end
    h = zeros(n,1);
    for i=1:niter
        h = h - tau * nablaF(h,K,y);
        Flist(i) = F(h,K,y);
    end
    % evaluate on a grid
    Theta = reshape( theta(kappa(G,X,sigma)*h) , [q,q]);
    % Display the classification probability.
    subplot(2,2,is);
    hold on; imagesc(t,t, Theta');
    options.ms = 5;
    plot_multiclasses(X,y,options); axis off; SetAR(1);
    colormap parula(256); caxis([0 1]);
    title(['\sigma=' num2str(sigma)]);
end
