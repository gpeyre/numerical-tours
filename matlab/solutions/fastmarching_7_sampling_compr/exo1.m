options.order = 2;
G = grad(M, options);
W = sum(G.^2,3);
% blur a little
W = perform_blurring(W, 10);
W = (W+epsilon).^alpha;
% scale to set up the contast
clf;
imageplot(M, 'Image', 1,2,1);
imageplot(W, 'Metric', 1,2,2);
