% centered kernel
J = eye(n) - ones(n)/n;
K = -1/2 * J*(D.^2)*J;
% diagonalization
opt.disp = 0; 
[Xstrain, val] = eigs(K, 2, 'LR', opt);
Xstrain = Xstrain .* repmat(sqrt(diag(val))', [n 1]);
Xstrain = Xstrain';
% plot graph
clf; hold on;
scatter(Xstrain(1,:),Xstrain(2,:),ms,v, 'filled'); 
plot_graph(A, Xstrain, options);
colormap jet(256);
axis('equal'); axis('off'); 
