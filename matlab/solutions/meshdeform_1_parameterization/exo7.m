f = 10;
v = cos(f*Y(1,:)*2*pi()) .* cos(f*Y(2,:)*2*pi());
% Display the function on the 2D parameteric domain.
options.face_vertex_color = rescale(v(:));
clf;
subplot(1,2,1);
colormap(jet(256));
plot_mesh([Y;zeros(1,n)],F,options);
colormap(jet(256));
view(2);
subplot(1,2,2);
clf;
colormap(jet(256));
plot_mesh(X,F,options);
colormap(jet(256));
