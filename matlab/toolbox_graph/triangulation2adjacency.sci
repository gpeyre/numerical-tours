function A = triangulation2adjacency(face,vertex)

// triangulation2adjacency - compute the adjacency matrix
//   of a given triangulation.
//
//   A = triangulation2adjacency(face);
//
//   Copyright (c) 2005 Gabriel Peyré


[tmp,face] = check_face_vertex([],face);
f = double(face)';

i = [f(:,1); f(:,1); f(:,2); f(:,2); f(:,3); f(:,3)];
j = [f(:,2); f(:,3); f(:,1); f(:,3); f(:,1); f(:,2)];
ij = [i,j];
v = ones(1,length(i))';
A = sparse(ij, v );
// avoid double links
// A = double(A>0);

endfunction