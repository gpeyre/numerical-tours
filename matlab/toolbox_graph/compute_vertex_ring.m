function vring = compute_vertex_ring(face)

% compute_vertex_ring - compute the 1 ring of each vertex in a triangulation.
%
%   vring = compute_vertex_ring(face);
%
%   vring{i} is the set of vertices that are adjacent
%   to vertex i.
%
%   Copyright (c) 2004 Gabriel Peyré

[tmp,face] = check_face_vertex([],face);

nverts = max(max(face));

A = triangulation2adjacency(face);
[i,j,s] = find(sparse(A));

% create empty cell array
vring{nverts} = [];

for m = 1:length(i)
    vring{i(m)}(end+1) = j(m);
end