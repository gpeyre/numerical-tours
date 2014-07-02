function [q,r] = compute_orthocenter(vertex,face)

% compute_orthocenter - compute the orthocenters of a triangulation
%
%   [q,r] = compute_orthocenter(vertex,face);
%
%   q are the centers.
%   r are the corresponding radii.
%
%   Copyright (c) 2008 Gabriel Peyre

% compute mid points
m1 = ( vertex(:,face(1,:)) + vertex(:,face(3,:)) )/2;
m2 = ( vertex(:,face(2,:)) + vertex(:,face(3,:)) )/2;
% compute normals
u1 = ( vertex(:,face(1,:)) - vertex(:,face(3,:)) );
u2 = ( vertex(:,face(2,:)) - vertex(:,face(3,:)) );
% compute determinant
d = u1(1,:).*u2(2,:) - u1(2,:).*u2(1,:);
d(abs(d)<1e-10)=1;
% compute inverse dual vectors <ui,vj>=delta(i-j)
v1 = [ u2(2,:)./d; -u1(2,:)./d];
v2 = [-u2(1,:)./d;  u1(1,:)./d];
% compute rhs
y = [ sum(m1.*u1); sum(m2.*u2)];
% solve system -> orthocenters
q = [ sum(v1.*y); sum(v2.*y) ];
% radius of inner circle
r = sqrt( sum( (q-vertex(:,face(1,:))).^2, 1 ) );