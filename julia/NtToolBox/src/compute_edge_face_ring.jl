function compute_edge_faces_ring(faces)
    """
        compute_edge_face_ring - compute faces adjacent to each edge

          e2f = compute_edge_face_ring(faces);

          e2f(i,j) and e2f(j,i) are the number of the two faces adjacent to
          edge (i,j).

          Copyright (c) 2007 Gabriel Peyre
    """
    n = maximum(faces)
    m = size(faces,2)

    i = [faces[1,:] faces[2,:] faces[3,:]]
    j = [faces[2,:] faces[3,:] faces[1,:]]
    s = [1:m 1:m 1:m];

    # first without duplicate
    tmp = unique(i + (maximum(i) + 1)*j) # unique
    I = indexin(tmp, i + (maximum(i) + 1)*j)

    # remaining items
    J = setdiff(1:length(s), I);

    # flip the duplicates
    i1 = [i[I]; j[J]]
    j1 = [j[I]; i[J]]
    s  = [s[I]; s[J]]

    # remove doublons
    tmp = unique(i1 + (maximum(i1) + 1)*j1)
    I = indexin(tmp, i1 + (maximum(i1) + 1)*j1)
    i1 = i1[I]
    j1 = j1[I]
    s  = s[I]
    A = sparse(i1,j1,s,n,n)

    # add missing points
    I = find(A'.!=0)
    I = I[A[I].==0]
    A[I]=-1

    return A
end;
