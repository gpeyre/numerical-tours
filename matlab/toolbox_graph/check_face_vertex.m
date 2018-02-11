function [vertex,face] = check_face_vertex(vertex,face, options)

% check_face_vertex - check that vertices and faces have the correct size
%
%   [vertex,face] = check_face_vertex(vertex,face);
%
%   Copyright (c) 2007 Gabriel Peyre

vertex = check_size(vertex,2,4);
face = check_size(face,3,4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a = check_size(a,vmin,vmax)
if isempty(a)
    return;
end
if size(a,1)>size(a,2)
    a = a';
end
if size(a,1)<3 && size(a,2)==3
    a = a';
end
if size(a,1)<=3 && size(a,2)>=3 && sum(abs(a(:,3)))==0
    % for flat triangles
    % a = a';
end
if size(a,1)<vmin ||  size(a,1)>vmax
    error('face or vertex is not of correct size');
end
