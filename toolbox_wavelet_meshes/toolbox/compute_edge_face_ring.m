function A = compute_edge_face_ring(face)

% compute_edge_face_ring - compute faces adjacent to each edge
%
%   e2f = compute_edge_face_ring(face);
%
%   e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to
%   edge (i,j).
%
%   Copyright (c) 2007 Gabriel Peyre


[tmp,face] = check_face_vertex([],face);

n = max(face(:));
m = size(face,2);
i = [face(1,:) face(2,:) face(3,:)];
j = [face(2,:) face(3,:) face(1,:)];
s = [1:m 1:m 1:m];

% first without duplicate
[tmp,I] = unique( i+(max(i)+1)*j );
% remaining items
J = setdiff(1:length(s), I);

% flip the duplicates
i1 = [i(I) j(J)];
j1 = [j(I) i(J)];
s = [s(I) s(J)];
A = sparse(i1,j1,s,n,n);


% add missing points
I = find( A'~=0 );
I = I( A(I)==0 ); 
A( I ) = -1;
