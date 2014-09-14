options.order = 2;
G = grad(f,options);
G = sqrt(sum(G.^2,3));
G = perform_blurring(G,3);
G = min(G,.4);
W = rescale(-G,.4,1);
clf;
imageplot(W);
