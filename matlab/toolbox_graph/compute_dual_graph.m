function [A,vertex1] = compute_dual_graph(face,vertex)

% compute_dual_graph - compute the dual graph of a given triangulation
%
%   [A,vertex1] = compute_dual_graph(face,vertex);
%
%   'A' is the adjacency matrix of the abstract dual graph
%   (recall that this graph link togeter adjacent faces
%   in the triangulation).
%
%   'vertex' is optional, and if given, the position of the vertex
%   of the dual graph (contained in 'vertex1') will 
%   the centroids of the face's vertex positions.
%   
%   Copyright (c) 2004 Gabriel Peyré

if nargin<2
    vertex = [];
end
[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,2);
nvert = max(max(face));

fring = compute_face_ring(face);

% compute the center of the faces 
if nargin>1
    vertex1 = [   ...
    sum(reshape(vertex(1,face),[3 nface]), 1)/3; ...
    sum(reshape(vertex(2,face),[3 nface]), 1)/3; ...
    sum(reshape(vertex(3,face),[3 nface]), 1)/3 ];
else
    vertex1 = [];
end

A = zeros(nface,nface);
for i=1:nface
    ring = fring{i};
    for j=1:length(ring)
        A(i,ring(j)) = 1;
        A(ring(j),i) = 1;
    end
end