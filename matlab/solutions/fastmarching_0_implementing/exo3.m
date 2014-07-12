options.method = 'fm';
D = perform_dijkstra_fm(W, x0, options);
k = 8;
displ = @(D)cos(2*pi*k*D/max(D(:)));
clf;
imageplot(displ(D));
colormap jet(256);
