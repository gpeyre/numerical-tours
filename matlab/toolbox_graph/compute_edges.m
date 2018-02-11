function edges = compute_edges(face)

% compute_edges - from a list of faces, compute the list of (unique) edges.
%
%   edge = compute_edges(face);
%
%   works for triangular and tets meshes.
%
%   Copyright (c) 2004 Gabriel Peyré

if isempty(face)
    edges=[];
    return;
end

[tmp,face] = check_face_vertex([],face);

if size(face,1)~=3 && size(face,1)~=4
    error('Problem, works for triangles and tets only.');
end

d = size(face,1);

edges = [];
for i=1:d
    sel = [i, mod(i,d)+1];
    edges = [edges, face(sel,:)];
end

% sort pair of vertex
I = find(edges(1,:)>edges(2,:));
J = find(edges(1,:)<=edges(2,:));
edges = [edges(end:-1:1,I), edges(:,J)];

% unique id
m = max(edges(:))+100;
id = edges(1,:) + m*edges(2,:);

[tmp,I] = unique(id);
edges = edges(:,I);