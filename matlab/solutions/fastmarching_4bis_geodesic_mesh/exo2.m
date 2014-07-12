options.niter = 200;
options.U = sqrt(sum(vertex.^2))';
[U,err] = perform_geodesic_iterative(vertex, faces, W, I, options);
clf;
options.face_vertex_color = mycolor(U, 3);
plot_mesh(vertex,faces,options);
colormap jet(256); shading interp;
