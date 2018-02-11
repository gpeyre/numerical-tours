function A = compute_triangulation_angles(vertex,face)

% compute_triangulation_angles - compute the vector of minimum angles
%
%   A = compute_triangulation_angles(vertex,face);
%
%   A(i) is the minium angle of face i.
%
%   Copyright (c) 2008 Gabriel Peyre

m = size(face,2);

A = zeros(3,m);

for i=1:3
    j1 = mod(i,3)+1;
    j2 = mod(i-2,3)+1;
    v1 = vertex(:,face(j1,:)) - vertex(:,face(i,:));
    v1 = v1 ./ repmat( sqrt(sum(v1.^2)), [2 1] );
    v2 = vertex(:,face(j2,:)) - vertex(:,face(i,:));
    v2 = v2 ./ repmat( sqrt(sum(v2.^2)), [2 1] );
    A(i,:) = acos( sum( v1.*v2, 1 ) );
end
A = min(A);