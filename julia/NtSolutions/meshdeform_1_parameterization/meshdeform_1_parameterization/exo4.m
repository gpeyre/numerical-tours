q = 64;
M = zeros(q,q,3);
for i=1:3
    M(:,:,i) = compute_triang_interp(F,Y,X(i,:), q);
end
[Y,X] = meshgrid(1:q,1:q);
T = mod(X+Y,2);
clf;
colormap(gray(256));
plot_surf_texture(M, T);
view(-40,70); zoom(1.5);
axis tight; axis square; axis off;
camlight;
