sigma_list = [5 6 8 10];
clf;
for is=1:length(sigma_list)
    sigma = sigma_list(is);
    kappa = @(X,Z)exp( -distmat(X,Z)/(2*sigma^2) );    
    %
    K = kappa(X,X);
    F = @(h)L(K*h,y);
    nablaF = @(h)K'*nablaL(K*h,y);
    % grad descent
    h = zeros(p,1);
    Flist = [];
    tau = 1000;
    niter = 200;
    for i=1:niter
        h = h - tau * nablaF(h);
        Flist(i) = F(h);
    end
    %
    Theta = theta(kappa(G,X)*h);
    Theta = reshape(Theta, [q,q]);
    % Display the classification probability.
    subplot(2,2,is);
    hold on;
    imagesc(t,t, Theta');
    for i=1:2
        I = find(y==2*i-3);
        plot(Z(I,1), Z(I,2), '.', 'Color', col(:,i), 'MarkerSize', ms);
    end
    colormap jet(256);
    caxis([0 1]);
    axis tight; axis equal; box on; axis off;
    title(['\sigma=' num2str(sigma)]);
end
