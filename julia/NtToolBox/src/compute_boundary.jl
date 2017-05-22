include("compute_edge_face_ring.jl")

function compute_boundary(face)
    """
        compute_boundary - compute the vertices on the boundary of a 3D mesh

          boundary=compute_boundary(face);

          Copyright (c) 2007 Gabriel Peyre
    """

    # compute edges (i,j) that are adjacent to only 1 face
    A = compute_edge_faces_ring(face)
    (i,j,v) = findnz(A)

    i = i[v .== -1]
    j = j[v .== -1]

    # build the boundary by traversing the edges
    boundary = [i[1]]
    i = i[2:end]
    j = j[2:end]
    while !(isempty(i))
        b = boundary[end]
        I = find(i.==b)
        if isempty(I)
            I = find(j.==b)
            if isempty(I)
                warn("Problem with boundary")
            end
            append!(boundary, i[I])
        else
            append!(boundary, j[I])
        end
        deleteat!(i, I)
        deleteat!(j, I)
    end

    #cyclic boundary
    #boundary = boundary + [boundary[0]]

    return boundary
end
