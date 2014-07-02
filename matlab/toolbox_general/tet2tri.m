function face = tet2tri(facet, vertex, keep_surface)

% tet2tri - convert a tet mesh to a tri mesh
%
%   face = tet2tri(facet, vertex, keep_surface);
%
%   if keep_surface==1, then keep only the outer part of the tet mesh.
%
%   Copyright (c) 2008 Gabriel Peyre

if nargin<3
    keep_surface = 0;
end

if size(facet,1)<size(facet,2)
    facet = facet';
end
if size(vertex,1)<size(vertex,2)
    vertex = vertex';
end

face=[facet(:,[1,2,3]);
    facet(:,[1,2,4]);
    facet(:,[1,3,4]);
    facet(:,[2,3,4])];
node4=[facet(:,4);facet(:,3);facet(:,2);facet(:,1)];

if keep_surface
    face=sort(face,2);
    [foo,ix,jx]=unique(face,'rows');
    vec=histc(jx,1:max(jx));
    qx=find(vec==1);
    face=face(ix(qx),:);
    node4=node4(ix(qx));
end

if not(isempty(vertex))
    % perform re-orientation
    v1=vertex(face(:,2),:)-vertex(face(:,1),:);
    v2=vertex(face(:,3),:)-vertex(face(:,1),:);
    v3=vertex(node4,:)-vertex(face(:,1),:);
    ix=find(dot(cross(v1,v2,2),v3,2)>0);
    face(ix,[2,3])=face(ix,[3,2]);
end

face = face';