def exo1():
    """
    Perform iterative smoothing and projection.
    Record the evolution of the number of inverted triangle in
    |ninvert|. Record also the evolution of the Dirichlet energy in
    |Edir|.
    """
    vertex1 = vertex
    vertex1 = vertex1 - repmat(mean(vertex1, 2), [1 n])
    vertex1 = vertex1 ./ repmat(sqrt(sum(vertex1.^2, 1)), [3 1])
    niter = 500
    ndisp = round([0 0.1 0.3 1]*niter); ndisp = max(ndisp, 1)
    k = 1
    ninvert = []
    Edir = []
    for i in 1: niter:
        % smooth 
        vertex1 = vertex1*tW'
        % project
        vertex1 = vertex1 ./ repmat(sqrt(sum(vertex1.^2, 1)), [3 1])
        % record Dirichlet energy
        E = 0
    for j in 1: 3:
            j1 = mod(j, 3) + 1
            % directed edge
            u = vertex1(: , faces(j, : )) - vertex1(: , faces(j1, : ))
            % norm squared
            u = sum(u.^2)
            % weights between the vertices
            w = W(faces(j, : ) + (faces(j1, : )-1)*n)
            E = sum(w.*u)
        Edir(end + 1) = E
        % record number of inverted triangles
        [normal, normalf] = compute_normal(vertex1, faces)
        C = squeeze(mean(reshape(vertex1(: , faces), [3 3 m]), 2))
        I = sum(C.*normalf)
        ninvert(end + 1) = sum(I <0)
        if i = =ndisp(k)
            % display
            subplot(2, 2, k)
            options.face_vertex_color = double(I(: ) >0)
            plot_mesh(vertex1, faces, options)
            colormap gray(256); axis tight
            shading faceted
            k = k + 1


def exo2():
    """
    Implement the mesh morphing.
    """


def exo3():
    """
    Perform the full descent.
    Record the decay of the energy in |Edir|.
    """
    vertex1 = vertex
    vertex1 = vertex1 - repmat(mean(vertex1, 2), [1 n])
    vertex1 = vertex1 ./ repmat(sqrt(sum(vertex1.^2, 1)), [3 1])
    niter = 500
    eta = .5
    ndisp = round([0 0.1 0.3 1]*niter); ndisp = max(ndisp, 1)
    k = 1
    Edir = []
    for it in 1: niter:
        % compute the center
        A = squeeze(sum(reshape(vertex1(: , faces), [3 3 m]), 2))
        % Compute the Dirichlet energy of each face.
        E = zeros(1, m)
    for i in 1: 3:
            i1 = mod(i, 3) + 1
            % directed edge
            u = vertex1(: , faces(i, : )) - vertex1(: , faces(i1, : ))
            % norm squared
            u = sum(u.^2)
            % weights between the vertices
            w = W(faces(i, : ) + (faces(i1, : )-1)*n)
            E = E + w.*u
        % Compute gradient direction.
        G = zeros(3, n)
        Edir(end + 1) = 0
    for j in 1: m:
            f = faces(: , j)
            Alpha = A(: , j)
            alpha = norm(Alpha)
    for i in 1: 3:
                i1 = mod(i  , 3) + 1
                i2 = mod(i + 1, 3) + 1
                % directed edges
                u1 = vertex1(: , f(i)) - vertex1(: , f(i1))
                u2 = vertex1(: , f(i)) - vertex1(: , f(i2))
                % weights between the vertices
                w1 = W(f(i) + (f(i1)-1)*n)
                w2 = W(f(i) + (f(i2)-1)*n)
                G(: , f(i)) = G(: , f(i)) + (w1*u1 + w2*u2) ./ alpha^2 - Alpha/ alpha^4 * E(j)
            Edir(end) = Edir(end) + E(j)/ alpha
        % Perform the gradient descent step and the projection.    
        vertex1 = vertex1 - eta*G
        vertex1 = vertex1 ./ repmat(sqrt(sum(vertex1.^2, 1)), [3 1])
        % display
        if 0
            if mod(it, 5) = =1
                plot_mesh(vertex1, faces, options)
                axis tight; shading faceted
                drawnow
        if i = =ndisp(k)
            % display
            subplot(2, 2, k)
            options.face_vertex_color = double(I(: ) >0)
            plot_mesh(vertex1, faces, options)
            colormap gray(256); axis tight
            shading faceted
            k = k + 1


