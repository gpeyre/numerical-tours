sigma_list = [.05 .1 .5 1 5 10];
%
clf;
for i=1:length(sigma_list)
    sigma = sigma_list(i);
    kappa = @(X,Z)exp( -distmat(X,Z)/(2*sigma^2) );
    % Regressor.
    h = (kappa(X,X)+lambda*eye(n))\y;
    Y = @(x)kappa(x,X)*h;
    % Eval on the grid
    yn = reshape(Y(Xn),[q,q]);
    %
    subplot(2,3,i);
    imagesc(t,t,yn); axis image; axis off; 
    colormap jet(256);
    title(['\sigma=' num2str(sigma)]);
end
