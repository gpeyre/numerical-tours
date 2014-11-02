def exo1():
    """
    Perform the same flattening, but with the combinatorial Laplacian.
    """
    L = compute_mesh_laplacian(vertex, faces, 'combinatorial', options)
    [U, S] = eig(full(L)); S = diag(S)
    [S, I] = sort(S, 'ascend'); U = U(: , I)
    vertexF = U(: , 2: 3)'
    vertexF = vertexF - repmat(vertexF(: , icenter), [1 n])
    theta = -pi/ 2 + atan2(vertexF(2, irotate), vertexF(1, irotate))
    vertexF = [vertexF(1, : )*cos(theta) + vertexF(2, : )*sin(theta); ...
               -vertexF(1, : )*sin(theta) + vertexF(2, : )*cos(theta)]
    plot_mesh(vertexF, faces)


def exo2():
    """
    Compute the embedding using Stress minimization with SMACOF.
    See the numerical tours on bending invariants for more details.
    """
    niter = 150
    stress = []
    vertexS = vertexF
    ndisp = [1 5 10 min(niter, 100) Inf]
    k = 1
    for i in 1: niter:
        if ndisp(k) = =i
            subplot(2, 2, k)
            plot_mesh(vertexS, faces, options)
            k = k + 1
        % Compute the distance matrix.
        D1 = repmat(sum(vertexS.^2, 1), n, 1)
        D1 = sqrt(D1 + D1' - 2*vertexS'*vertexS)
        % Compute the scaling matrix.
        B = -D./ max(D1, 1e-10)
        B = B - diag(sum(B))
        % update
        vertexS = (B*vertexS')' / n
        % Xstress = Xstress-repmat(mean(Xstress, 2), [1 n])
        % record stress
        stress(end + 1) = sqrt(sum(abs(D(: )-D1(: )).^2) / n^2)


def exo3():
    """
    Compute mesh parameterization using a circle as boundary.
    otan weights
    oundary
    ixed positions
    ystem
    et up the right hand sizes with the fixed position.
    olve
    lign
    isplay
    """
    options.symmetrize = 1
    options.normalize = 0
    L = compute_mesh_laplacian(vertex, faces, 'conformal', options)
    options.verb = 0
    boundary = compute_boundary(faces, options)
    p = length(boundary)
    t = linspace(0, 2*pi, p + 1)'; t(p) = []
    x0 = cos(t); y0 = sin(t)
    L1 = L
    L1(boundary, : ) = 0
    L1(boundary + (boundary-1)*n) = 1
    Rx = zeros(n, 1); Rx(boundary) = x0
    Ry = zeros(n, 1); Ry(boundary) = y0
    x = L1 \ Rx
    y = L1 \ Ry
    vertexF = [x'; y']
    vertexF = vertexF - repmat(vertexF(: , icenter), [1 n])
    theta = -pi/ 2 + atan2(vertexF(2, irotate), vertexF(1, irotate))
    vertexF = [vertexF(1, : )*cos(theta) + vertexF(2, : )*sin(theta); ...
               -vertexF(1, : )*sin(theta) + vertexF(2, : )*cos(theta)]
    plot_mesh(vertexF, faces, options)


