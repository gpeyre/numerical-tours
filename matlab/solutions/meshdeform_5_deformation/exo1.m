Delta0 = zeros(3,n);
Delta0(3,I) = vertex(1,I).^2 - vertex(2,I).^2;
% deform
vertex1 = vertex+( L1 \ Delta0' )';
% Display it.
clf;
plot_mesh(vertex1,faces);
view(-150,45);
