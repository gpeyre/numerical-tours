W = rescale( min(C,.1), .001, 1);
options.W = W;
%
landmarks = [5000];
options.constraint_map = [];
[D,Z,Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options);
% display
clf;
hold on;
options.face_vertex_color = mod( 20*D/max(D),1 );
plot_mesh(vertex,faces, options);
colormap jet(256);
h = plot3(vertex(1,landmarks), vertex(2,landmarks), vertex(3,landmarks), 'r.');
set(h, 'MarkerSize', 20);
