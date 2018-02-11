function [normal,normalf] = compute_normal(vertex,face)

% compute_normal - compute the normal of a triangulation
%
%   [normal,normalf] = compute_normal(vertex,face);
%
%   normal(i,:) is the normal at vertex i.
%   normalf(j,:) is the normal at face j.
%
%   Copyright (c) 2004 Gabriel Peyré

[vertex,face] = check_face_vertex(vertex,face);

nface = size(face,1);
nvert = size(vertex,1);
normal = zeros(nvert,3);

% unit normals to the faces
normalf = crossp( vertex(face(:,2),:)-vertex(face(:,1),:), ...
                  vertex(face(:,3),:)-vertex(face(:,1),:) );
d = sqrt( sum(normalf.^2,2) ); d(d<eps)=1;
normalf = normalf ./ repmat( d, 1,3 );

% unit normal to the vertex
normal = zeros(nvert,3);
for i=1:nface
    f = face(i,:);
    for j=1:3
        normal( f(j),: ) = normal( f(j),: ) + normalf(i,:);
    end
end
% normalize
d = sqrt( sum(normal.^2,2) ); d(d<eps)=1;
normal = normal ./ repmat( d, 1,3 );

% enforce that the normal are outward
v = vertex - repmat(mean(vertex,2), 1,3);
s = sum( v.*normal, 1 );
if sum(s>0)<sum(s<0)
    % flip
    normal = -normal;
    normalf = -normalf;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function z = crossp(x,y)
% x and y are (m,3) dimensional
z = x;
z(:,1) = x(:,2).*y(:,3) - x(:,3).*y(:,2);
z(:,2) = x(:,3).*y(:,1) - x(:,1).*y(:,3);
z(:,3) = x(:,1).*y(:,2) - x(:,2).*y(:,1);
