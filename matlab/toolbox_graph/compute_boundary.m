function boundary=compute_boundary(face, options)

% compute_boundary - compute the vertices on the boundary of a 3D mesh
%
%   boundary=compute_boundary(face);
%
%   Copyright (c) 2007 Gabriel Peyre

if size(face,1)<size(face,2)
    face=face';
end

%% compute edges (i,j) that are adjacent to only 1 face
A = compute_edge_face_ring(face);
[i,j,v] = find(A);
i = i(v==-1);
j = j(v==-1);

%% build the boundary by traversing the edges
boundary = i(1); i(1) = []; j(1) = [];
while not(isempty(i))
    b = boundary(end);
    I = find(i==b);
    if isempty(I)
        I = find(j==b);
        if isempty(I)
            warning('Problem with boundary');
            break;
        end
        boundary(end+1) = i(I);
    else
        boundary(end+1) = j(I);
    end        
	i(I) = []; j(I) = [];
end

return;

%% OLD CODE %%

nvert=max(max(face));
nface=size(face,1);


% count number of faces adjacent to a vertex
A=sparse(nvert,nvert);
for i=1:nface
    if verb
        progressbar(i,nface);
    end
    f=face(i,:);
    A(f(1),f(2))=A(f(1),f(2))+1;
    A(f(1),f(3))=A(f(1),f(3))+1;
    A(f(3),f(2))=A(f(3),f(2))+1;
end
A=A+A';

for i=1:nvert
    u=find(A(i,:)==1);
    if ~isempty(u)
        boundary=[i u(1)];
        break;
    end
end

s=boundary(2);
i=2;
while(i<=nvert)
    u=find(A(s,:)==1);
    if length(u)~=2
        warning('problem in boundary');
    end
    if u(1)==boundary(i-1)
        s=u(2);
    else
        s=u(1);
    end
    if s~=boundary(1)
        boundary=[boundary s];
    else
        break;
    end
    i=i+1;
end
       
if i>nvert
    warning('problem in boundary');
end


%%% OLD %%%
function v = compute_boundary_old(faces)

nvert = max(face(:));
ring = compute_vertex_ring( face );

% compute boundary
v = -1;
for i=1:nvert   % first find a starting vertex
    f = ring{i};
    if f(end)<0
        v = i;
        break;
    end
end
if v<0
    error('No boundary found.');
end
boundary = [v];
prev = -1;
while true
    f = ring{v};
    if f(end)>=0
        error('Problem in boundary');
    end
    if f(1)~=prev
        prev = v;
        v = f(1);
    else
        prev = v;
        v = f(end-1);
    end
    if ~isempty( find(boundary==v) )
        % we have reach the begining of the boundary
        if v~=boundary(1)
            warning('Begining and end of boundary doesn''t match.');
        else
            break;
        end
    end
    boundary = [boundary,v];
end