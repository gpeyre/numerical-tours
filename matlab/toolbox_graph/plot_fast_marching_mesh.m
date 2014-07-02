function plot_fast_marching_mesh(vertex,faces, D, paths, options)

% plot_fast_marching_mesh - plot the result of the fast marching on a mesh.
%
%   plot_fast_marching_mesh(vertex,faces, D, paths, options);
%
%   paths can be a cell array of path (to display multiple geodesics).
%   paths should be computed with compute_geodesic_mesh.
%
%   options.voronoi_edges is a list of voronoi edges computed with
%       compute_voronoi_mesh.
%
%   See also perform_fast_marching_mesh, compute_geodesic_mesh, compute_voronoi_mesh.
%
%   You can force color equalization using options.colorfx = 'equalize' 
%   options.start_points is a set of start points (in red), for instance the center of geodesic cells.
%
%   Copyright (c) 2007 Gabriel Peyre

options.null = 0;
if not(iscell(paths))
    paths = {paths};
end

col = D;
col(col==Inf) = 0;

start_points = getoptions(options, 'start_points', []);
ve = getoptions(options, 'voronoi_edges', []);
ms = getoptions(options, 'ms', 20);
lw = getoptions(options, 'lw', 3);

colorfx = getoptions(options, 'colorfx', 'none');
if strcmp(colorfx, 'wrapping')
    num = 20;
    col = rescale(col);
    col = mod(col*num,2); col(col>1) = 2-col(col>1);
end
if strcmp(colorfx, 'equalize')    
    col(col>0) = perform_hist_eq(col(col>0), 'linear');
end

options.face_vertex_color = col;
% clf;
hold on;
plot_mesh(vertex, faces, options);
for i=1:length(paths)
    path = paths{i};
    if not(isempty(path))
        end_point = path(:,1);
        start_point = path(:,end);
        % display starting points
        h = plot3( start_point(1),start_point(2), start_point(3), 'r.');
        set(h, 'MarkerSize', ms);
        % display endding points
        h = plot3( end_point(1), end_point(2), end_point(3), 'b.');
        set(h, 'MarkerSize', ms);
        % display geodesic
        h = plot3(path(1,:), path(2,:),path(3,:), 'k');
        set(h,'LineWidth', lw);
    end
end
% display starting points
for i=1:length(start_points)
    start_point = vertex(:,start_points(i));
    h = plot3( start_point(1),start_point(2), start_point(3), 'r.');
    set(h, 'MarkerSize', ms);    
end
% display voronoi
for i=1:size(ve,2)
    h = plot3( [ve(1,i) ve(4,i)], ...
        [ve(2,i) ve(5,i)], ...
        [ve(3,i) ve(6,i)], 'r');
    set(h, 'LineWidth', lw);
end
hold off;
cm = jet(256);
cm(1,:) = [1 1 1]/5;
colormap(cm);
% camlight;