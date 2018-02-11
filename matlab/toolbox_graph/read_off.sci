function [vertex,face] = read_off(filename)

// read_off - read data from OFF file.
//
//   [vertex,face] = read_off(filename);
//
//   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices.
//   'face' is a 'nb.face x 3' array specifying the connectivity of the mesh.
//
//   Copyright (c) 2003 Gabriel Peyré


fid = mopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

// str = fgets(fid);   // -1 if eof
[cnt,str] = mfscanf(1, fid,'%s');
if 0 // ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

[cnt,nvert,nface,u] = mfscanf(1, fid, '%d %d %d');

// read vertices
[cnt,vertex] = mfscanf(nvert*3, fid,'%f');
if length(vertex)~=3*nvert
    warning('Problem in reading vertices.');
end
vertex = matrix(vertex,3,nvert);

// read Face 1  1088 480 1022
[cnt, face] = mfscanf(nface*4, fid, '%d');
if length(face)~=4*nface
    warning('Problem in reading faces.');
end
face = matrix(face, 4, nface);
face = face(2:4,:)+1;

mclose(fid);

endfunction