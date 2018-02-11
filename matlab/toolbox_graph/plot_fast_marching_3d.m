function plot_fast_marching_3d(W,S,path,start_points,end_points, options)

% plot_fast_marching_3d - plot the result of the fast marching.
%
%   plot_fast_marching_3d(W,S,path [,start_points,end_points, options] );
%
%   If you provide W, path is assumed to lie in [1,n]^3 where n=length(W).
%   Path can be a cell array of 3D curves or a single 3D curve.
%
%   Copyright (c) 2004 Gabriel Peyré

options.null = 0;
if isfield(options, 'plot_isosurface')
    plot_isosurface = options.plot_isosurface;
else
    plot_isosurface = 1;
end  
if isfield(options, 'plot_planes')
    plot_planes = options.plot_planes;
else
    plot_planes = 1;
end  
if isfield(options, 'path_width')
    path_width = options.path_width;
else
    path_width = 5;
end
if isfield(options, 'point_size')
    point_size = options.point_size;
else
    point_size = 16;
end


[n,p,q] = size(W);

hold on;

if ~isempty(W) && plot_planes
    if 0
    open_list = find( S==0 );
    close_list = find( S==-1 );
    W(open_list) = 0;
    W(close_list) = W(close_list) - 0.2;
    end
    plot_volumetric_data(W);
end

if nargin>=4 && ~isempty(start_points)
    if size(start_points,1)~=3
        start_points = start_points';
    end
    if size(start_points,1)~=3
        error('start_points must be of size 3xk');
    end
%    start_points = (start_points-1)/(n-1);
    start_points = shift_array(start_points);
    h = plot3(start_points(1), start_points(2), start_points(3), 'o', 'MarkerSize', point_size, ...
                    'MarkerEdgeColor','k', 'MarkerFaceColor','b');
end

if nargin>=5 && ~isempty(end_points)
    if size(end_points,1)~=3
        start_points = start_points';
    end
    if size(end_points,1)~=3
        error('start_points must be of size 2xk');
    end
    end_points = shift_array(end_points);
%    end_points = (end_points-1)/(n-1);
    plot3(end_points(1), end_points(2), end_points(3), 'o', 'MarkerSize', point_size, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','r');
end

if nargin>=3 && ~isempty(path)
    if ~iscell(path)
        path = {path};
    end
    for i=1:length(path)
		path1 = path{i};
        if size(path1,2)~=3
            path1 = path1';
        end
        path1 = shift_array(path1);
        plot3(path1(:,1), path1(:,2), path1(:,3), '-', 'LineWidth', path_width);
    end
end

if plot_isosurface
    S = smooth3(S,'gaussian',5,1);
    F = isosurface( S,0.5 );
    p = patch(F);
    isonormals( S,p );
    set(p, 'FaceColor', 'red', 'EdgeColor', 'none'); 
end

box on;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'ZTick', []);
lighting phong;
alpha(0.7);    
camproj('perspective');
view(3);
axis equal;
axis([1 size(W,1) 1 size(W,2) 1 size(W,3)]);
daspect([1 1 1]);
% cameramenu;
camlight;

hold off;

function path1 = shift_array(path)
if( size(path, 2)==1 )
    path1 = path( [2:-1:1 end] );
else
    path1 = path( :, [2:-1:1 end] );
end

function plot_volumetric_data(W, nb_colors)

% plot_volumetric_data - plot a cube of data.
%
%   plot_volumetric_data(W, nb_colors);
%
%   Copyright (c) 2005 Gabriel Peyré

if nargin<2
    nb_colors = 256;
end

[n,p,q] = size(W);
h = slice(W,p/2,n/2,q/2);
set(h,'FaceColor','interp','EdgeColor','none');

box on;
set(gca, 'XTick', []);
set(gca, 'YTick', []);
set(gca, 'ZTick', []);

colormap( jet(nb_colors) );

lighting phong;
camlight infinite; 
camproj('perspective');

view(3);
axis tight;
axis equal;
% cameramenu;