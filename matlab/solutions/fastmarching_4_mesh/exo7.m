% compute distances
options.nb_iter_max = Inf;
options.W = ones(nvert,1);
[D0,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
options.W = W;
[D,S,Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options);
% compute paths
pend = [9969; 5073];
options.method = 'continuous';
paths0 = compute_geodesic_mesh(D0, vertex, faces, pend, options);
paths = compute_geodesic_mesh(D, vertex, faces, pend, options);
% display
clf;
subplot(1,2,1);
plot_fast_marching_mesh(vertex,faces, f, paths0, options); 
subplot(1,2,2);
plot_fast_marching_mesh(vertex,faces, f, paths, options); 
colormap gray(256);
%EXO