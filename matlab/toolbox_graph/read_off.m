function [vertex,face] = read_off(filename)

% read_off - read data from OFF file.
%
%   [vertex,face] = read_off(filename);
%
%   'vertex' is a 'nb.vert x 3' array specifying the position of the vertices in 3-D.
%   'face' is a 'nb.face x k' array specifying the connectivity of the mesh
%   where k=3 for triangular mesh, and k=4 for quadmesh.
%
%   Copyright (c) 2012 Gabriel Peyré


fid = fopen(filename,'r');
if( fid==-1 )
    error('Can''t open the file.');
    return;
end

str = fgets(fid);   % -1 if eof
if ~strcmp(str(1:3), 'OFF')
    error('The file is not a valid OFF one.');    
end

str = fgets(fid);
[a,str] = strtok(str); nvert = str2num(a);
[a,str] = strtok(str); nface = str2num(a);

[A,cnt] = fscanf(fid,'%f %f %f', 3*nvert);
if cnt~=3*nvert
    warning('Problem in reading vertices.');
end
A = reshape(A, 3, cnt/3);
vertex = A;
% read Face 1 1088 480 1022  (triangular mesh)
[k,cnt] = fscanf(fid,'%d', 1); % read wether triangular mesh or quad mesh
if k==3
    %% triangular mesh
    [a,cnt] = fscanf(fid,'%d %d %d\n', 3);
elseif k==4
    %% quad mesh
    [a,cnt] = fscanf(fid,'%d %d %d %d\n', 4);    
else
    error('Can only read triangular or quad meshes');    
end
[A,cnt] = fscanf(fid,'%d %d %d %d\n', (k+1)*(nface-1));
if cnt~=(k+1)*(nface-1)
    warning('Problem in reading faces.');
end
A = reshape(A, k+1, cnt/(k+1));
face = [a(:)+1, A(2:(k+1),:)+1];


fclose(fid);

