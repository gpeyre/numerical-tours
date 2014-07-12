clf;
h = vol3d('cdata',M,'texture','2D');
view(3); axis off;
% set up a colormap
colormap bone(256);
% set up an alpha map
options.sigma = .08; % control the width of the non-transparent region
options.center = .6; % here a value in [0,1]
a = compute_alpha_map('gaussian', options); % you can plot(a) to see the alphamap
% refresh the rendering
vol3d(h);
