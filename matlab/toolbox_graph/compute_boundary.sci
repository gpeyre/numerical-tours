function boundary=compute_boundary(faces, options)

// compute_boundary - compute the vertices on the boundary of a 3D mesh
//
//   boundary=compute_boundary(faces);
//
//   Copyright (c) 2007 Gabriel Peyre

if size(faces,1)<size(faces,2)
    faces=faces';
end


//// compute edges (i,j) that are adjacent to only 1 face
A = compute_edge_face_ring(faces);
// [i,j,v] = find(A);
// i = i(v==-1);
// j = j(v==-1);
[ij,v,mn]=spget(A);
i = ij(v==-1,1);
j = ij(v==-1,2);

//// build the boundary by traversing the edges
boundary = i(1); i(1) = []; j(1) = [];
while not(isempty(i))
    b = boundary($);
    I = find(i==b);
    if isempty(I)
        I = find(j==b);
        if isempty(I)
            warning('Problem with boundary');
            break;
        end
        boundary($+1) = i(I);
    else
        boundary($+1) = j(I);
    end        
	i(I) = []; j(I) = [];
end


endfunction

