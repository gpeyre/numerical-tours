q = 12;
vertex(1,2:q+1) = .97;
vertex(2,2:q+1) = linspace(0.03,.97,q);
faces = compute_delaunay(vertex);
W = ones(n,1);  W(vertex(1,:)<.5) = 1/2;
I = 1;
options.niter = 300;
options.U = sqrt( sum( (vertex-repmat(vertex(:,I),[1 n])).^2 ) )';
[U,err] = perform_geodesic_iterative(vertex, faces, W, I, options);
%
options.method = 'continuous';
paths = compute_geodesic_mesh(U, vertex, faces, 2:q+1, options); 
clf;
plot_fast_marching_mesh(vertex,faces, mycolor(U, 8), paths, options);
%EXO