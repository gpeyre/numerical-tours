function plot_triangulation(vertex,faces, M, options)

% plot_triangulation - plot a 2D triangulation
%
%   plot_triangulation(vertex,face,M, options);
%
%   The point are assume to be in (1,n).
%
%   Copyright (c) 2009 Gabriel Peyre

options.null = 0;
lw = getoptions(options, 'edgewidth', 2);
ms = getoptions(options, 'vertexsize', 20);

edges = compute_edges(faces);

vertex = vertex(2:-1:1,:);

if not(isempty(M))
    imageplot(M);
end
hold on;
hh = plot_edges(edges, vertex, 'k');
set(hh, 'LineWidth', lw);
hh = plot(vertex(1,:),vertex(2,:), 'r.');
set(hh, 'MarkerSize', ms);
hold off;

function h = plot_edges(edges, vertex, color)

% plot_edges - plot a list of edges
%
%   h = plot_edges(edges, vertex, color);
%
%   Copyright (c) 2004 Gabriel Peyré

if nargin<3
    color = 'b';
end

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end

x = [ vertex(1,edges(1,:)); vertex(1,edges(2,:)) ];
y = [ vertex(2,edges(1,:)); vertex(2,edges(2,:)) ];
if size(vertex,1)==2
    h = line(x,y, 'color', color);
elseif size(vertex,1)==3
    z = [ vertex(3,edges(1,:)); vertex(3,edges(2,:)) ];
    h = line(x,y,z, 'color', color);    
else
    error('Works only for 2D and 3D plots');    
end