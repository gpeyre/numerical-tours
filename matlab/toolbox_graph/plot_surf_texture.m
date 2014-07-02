function h = plot_surf_texture(M, T)

% plot_surf_texture - plot a surface with a texture on it.
%
%   h = plot_surf_texture(M, T);
%
%   M(:,:,i) for i=1,2,3 are the 3D coordonate of the surface.
%   T is a 2D image mapped on the surface.
%
%   Copyright (c) 2010 Gabriel Peyre

h = surf(M(:,:,1), M(:,:,2), M(:,:,3) );
set(h,'CData',rescale(T')*255,'FaceColor','texturemap', 'EdgeColor','none', 'CDataMapping','direct');
colormap(gray(256));
