function ring = compute_vertex_face_ring(face)

% compute_vertex_face_ring - compute the faces adjacent to each vertex
%
%   ring = compute_vertex_face_ring(face);
%
%   Copyright (c) 2007 Gabriel Peyr?


[tmp,face] = check_face_vertex([],face);

nfaces = size(face,2);
nverts = max(face(:));

ring{nverts} = [];

for i=1:nfaces
    for k=1:3
        ring{face(k,i)}(end+1) = i;
    end
end