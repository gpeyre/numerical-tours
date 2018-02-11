function face = compute_delaunay(vertex)

% compute_delaunay - compute faces of a Delaunay triangulation
%
%   face = compute_delaunay(vertex);
%
%   Copyright (c) 2008 Gabriel Peyre

if size(vertex,1)>size(vertex,2)
    vertex = vertex';
end

vertex = vertex + randn(size(vertex))*1e-9;

warning off;
face = delaunay(vertex(1,:), vertex(2,:))';
warning on;