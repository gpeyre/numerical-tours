def exo1():
    """
    Perform the full subdivision.
    """
    face = face0
    vertex = vertex0
    for i in 1: 4:
        edge = compute_edges(face)
        n = size(vertex, 2)
        ne = size(edge, 2)
        % Compute the number of the three edges associated to each face.
        A = sparse([edge(1, : ); edge(2, : )], [edge(2, : ); edge(1, : )], [n + (1: ne); n + (1: ne)], n, n)
        v12 = full(A(face(1, : ) + (face(2, : )-1)*n))
        v23 = full(A(face(2, : ) + (face(3, : )-1)*n))
        v31 = full(A(face(3, : ) + (face(1, : )-1)*n))
        % Compute the new faces, each old face generates 4 faces.
        face = [cat(1, face(1, : ), v12, v31), ...
            cat(1, face(2, : ), v23, v12), ...
            cat(1, face(3, : ), v31, v23), ...
            cat(1, v12, v23, v31)]
        % Add new vertices at the edges center.
        vertex = [vertex, (vertex(: , edge(1, : )) + vertex(: , edge(2, : )))/ 2]
        % Project the points on the sphere.
        d = sqrt(sum(vertex.^2, 1))
        vertex = vertex ./ repmat(d, [size(vertex, 1) 1])
        % display
        subplot(2, 2, i)
        plot_mesh(vertex, face)
        shading('faceted'); lighting('flat'); view(3); axis('tight')


def exo2():
    """
    Try with other control meshes.
    """


def exo3():
    """
    Perform several steps of subdivision.
    """
    face = face0
    vertex = vertex0
    for i in 1: 2:
        % Perform the subdivision.
        edge = compute_edges(face)
        n = size(vertex, 2)
        ne = size(edge, 2)
        % Compute the number of the three edges associated to each face.
        A = sparse([edge(1, : ); edge(2, : )], [edge(2, : ); edge(1, : )], [n + (1: ne); n + (1: ne)], n, n)
        v12 = full(A(face(1, : ) + (face(2, : )-1)*n))
        v23 = full(A(face(2, : ) + (face(3, : )-1)*n))
        v31 = full(A(face(3, : ) + (face(1, : )-1)*n))
        % Compute the new faces, each old face generates 4 faces.
        face_old = face
        face = [cat(1, face(1, : ), v12, v31), ...
            cat(1, face(2, : ), v23, v12), ...
            cat(1, face(3, : ), v31, v23), ...
            cat(1, v12, v23, v31)]
        % Compute the vertex and face ring.
        global vring e2f fring facej
        vring = compute_vertex_ring(face)
        e2f = compute_edge_face_ring(face_old)
        fring = compute_face_ring(face_old)
        facej = face_old
        % Compute the interpolated position using
    for k in n + 1: n + ne:
            [e, v, g] = compute_butterfly_neighbors(k, n)
            vertex(: , k) = 1/ 2*sum(vertex(: , e), 2) + 1/ 8*sum(vertex(: , v), 2) - 1/ 16*sum(vertex(: , g), 2)
    subplot(1, 2, 1)
    plot_mesh(vertex0, face0, options)
    shading('faceted'); lighting('flat'); axis('tight')
    subplot(1, 2, 2)
    plot_mesh(vertex, face, options)
    shading('interp'); lighting('phong'); axis('tight')


def exo4():
    """
    Try on different 3D models.
    """


def exo5():
    """
    Implement another subdivision scheme that is not interpolating, for
    instance the loop scheme. Be careful about the handling of points that
    does not have valence 6.
    """
    face = face0
    vertex = vertex0
    for i in 1: 2:
        % Perform the subdivision.
        edge = compute_edges(face)
        n = size(vertex, 2)
        ne = size(edge, 2)
        % Compute the number of the three edges associated to each face.
        A = sparse([edge(1, : ); edge(2, : )], [edge(2, : ); edge(1, : )], [n + (1: ne); n + (1: ne)], n, n)
        v12 = full(A(face(1, : ) + (face(2, : )-1)*n))
        v23 = full(A(face(2, : ) + (face(3, : )-1)*n))
        v31 = full(A(face(3, : ) + (face(1, : )-1)*n))
        % Compute the new faces, each old face generates 4 faces.
        face_old = face
        face = [cat(1, face(1, : ), v12, v31), ...
            cat(1, face(2, : ), v23, v12), ...
            cat(1, face(3, : ), v31, v23), ...
            cat(1, v12, v23, v31)]
        % Compute the vertex and face ring.
        global vring e2f fring facej
        vring = compute_vertex_ring(face)
        vring0 = compute_vertex_ring(face_old)
        e2f = compute_edge_face_ring(face_old)
        fring = compute_face_ring(face_old)
        facej = face_old
        % move old vertices
        vertex1 = vertex
    for k in 1: n:
            m = length(vring0{k})
            beta = 1/ m*(5/ 8 - (3/ 8 + 1/ 4*cos(2*pi/ m))^2);   % loop original construction
            %    beta = 3/ (8*m);         % warren weights
            vertex1(: , k) = vertex(: , k)*(1-m*beta) + beta*sum(vertex(: , vring0{k}), 2)
        vertex = vertex1
        % move new vertices
    for k in n + 1: n + ne:
            [e, v] = compute_butterfly_neighbors(k, n)
            vertex(: , k) = 3/ 8*sum(vertex(: , e), 2) + 1/ 8*sum(vertex(: , v), 2)
    subplot(1, 2, 1)
    plot_mesh(vertex0, face0, options)
    shading('faceted'); lighting('flat'); axis('tight')
    subplot(1, 2, 2)
    plot_mesh(vertex, face, options)
    shading('interp'); lighting('phong'); axis('tight')


def exo6():
    """
    Implement another subdivision scheme that does not perform a 1:4 split
    of each face, for instance the sqrt(3) scheme.
    """


