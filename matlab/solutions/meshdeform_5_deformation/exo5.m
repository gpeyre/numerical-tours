normal1 = compute_normal(vertex1,faces);
vertex2 = vertex1 + d .* normal1;
% Display it.
clf;
plot_mesh(vertex2,faces);
view(-150,45);
