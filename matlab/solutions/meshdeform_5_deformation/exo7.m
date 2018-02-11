LL1 = L0*L0;
LL1(I,:) = 0;
LL1(I + (I-1)*n) = 1;
% deform
vertex1 = vertex0+( LL1 \ Delta0' )';
% Display it.
clf;
plot_mesh(vertex1,faces);
view(-150,45);
zoom(.8);
