function ring = compute_vertex_face_ring(face)

// compute_vertex_face_ring - compute the faces adjacent to each vertex
//
//   ring = compute_vertex_face_ring(face);
//
//   Copyright (c) 2007 Gabriel Peyre


[tmp,face] = check_face_vertex([],face);

nfaces = size(face,2);
nverts = max(face(:));

ring = cell(nverts,1);
ring(nverts).entries = [];

for i=1:nfaces
    for k=1:3
        a = ring(face(k,i)).entries;
		a($+1) = i;
		ring(face(k,i)).entries = a;
    end
end

endfunction