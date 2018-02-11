function v = compute_mesh_volume(vertex,face)

% compute_mesh_volume - compute the volume of an oriented mesh.
%
%   v = compute_mesh_volume(vertex,face);
%
%   Copyright (c) 2008 Gabriel Peyre

m = size(face,2);
A = zeros(3,3,m);
for i=1:3
    A(:,i,:) = reshape(vertex(:,face(i,:)), [3 1 m]);
end
v = abs(sum(det3(A))/4);