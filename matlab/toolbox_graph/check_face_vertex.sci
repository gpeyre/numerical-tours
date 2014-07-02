function [vertex,face] = check_face_vertex(vertex,face, options)

// check_face_vertex - check that vertices and faces have the correct size
//
//   [vertex,face] = check_face_vertex(vertex,face);
//
//   Copyright (c) 2007 Gabriel Peyre

vertex = check_size(vertex);
face = check_size(face);

endfunction

//////////////////////////////////////////////////////////
function a = check_size(a)

if isempty(a)
    return;
end
if size(a,1)>size(a,2)
    a = a';
end
if size(a,1)<3 & size(a,2)==3
    a = a';
end
if size(a,1)<=3 & size(a,2)>=3 & sum(abs(a(:,3)))==0
    // for flat triangles
    a = a';
end
if size(a,1)~=3 & size(a,1)~=4
    error('face or vertex is not of correct size');
end

endfunction