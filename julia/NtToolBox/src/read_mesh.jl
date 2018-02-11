function read_mesh(name)
    """
        reading from a OFF file in 3 dimensions, returning X0 (coordinates) and F (faces)
    """
    file = open(name)

    #check type of file
    file_type = strip(readline(file))
    if file_type != "OFF"
        throw(ArgumentError("Wrong type of file, only reads OFF files"))
    end

    #number of vertices/faces/edges:
    (n_verts, n_faces, n_edges) = tuple([parse(Int,s) for s in split(strip(readline(file)),' ')]...)

    #vertices
    X0 = zeros(3,n_verts)
    for i in 1:n_verts
        X0[:,i]= [parse(Float64,s) for s in split(strip(readline(file)),' ')]
    end

    #faces
    F = zeros(Int,3,n_faces)
    for i in 1:n_faces
        F[:,i] = [parse(Int,s) for s in split(strip(readline(file)),' ')[2:end]]
    end

    return X0, F+1
end
