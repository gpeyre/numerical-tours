def exo1():
    """
    Compute a density function that is larger at area of large gradient.
    |W(x) = (norm(grad(M))+epsilon|)^alpha|, for |alpha=.7|.
    To stabilize the process, you can smooth a bit the gradient magnitude.
    lur a little
    cale to set up the contast
    """
    options.order = 2
    G = grad(M, options)
    W = sum(G.^2, 3)
    W = perform_blurring(W, 10)
    W = (W + epsilon).^alpha
    imageplot(M, 'Image', 1, 2, 1)
    imageplot(W, 'Metric', 1, 2, 2)


def exo2():
    """
    Perform farthest points sampling to compute sampling location |vertex|
    and the corresponding geodesic Delaunay triangulation |faces|.
    ompute the Delaunay triangulation
    """
    ndisp = [20 50 100 m]
    vertex = [1; 1]
    [D, Z, Q] = perform_fast_marching(1./ W, vertex)
    vertex_svg = {}
    faces_svg = {}
    clf; u = 1
    for k in 2: m:
        options.constraint_map = D
        [D1, Z, Q] = perform_fast_marching(1./ W, vertex(: , end), options)
        D = min(D, D1)
        [tmp, i] = max(D(: ))
        [x, y] = ind2sub([n n], i)
        vertex(: , end + 1) = [x; y]
        if k = =ndisp(u)
            subplot(2, 2, u)
            hold on
            imageplot(M, [num2str(k) ' points'])
            plot(vertex(2, : ), vertex(1, : ), 'r.')
            axis ij
            u = u + 1
    [D, Z, Q] = perform_fast_marching(1./ W, vertex)
    faces = compute_voronoi_triangulation(Q, vertex)


def exo3():
    """
    For a large value of |m| compute the approximation for several |alpha|.
    """
    alpha_list = [.4 .7 1 1.3]
    m = 600
    for ialpha in 1: length(alpha_list):
        alpha = alpha_list(ialpha)
        % metric
        W = sum(G.^2, 3)
        W = perform_blurring(W, 6)
        W = (W + epsilon).^alpha
        % Farthest point
        vertex = [1; 1]
        [D, Z, Q] = perform_fast_marching(1./ W, vertex)
    for k in 2: m:
            options.constraint_map = D
            [D1, Z, Q] = perform_fast_marching(1./ W, vertex(: , end), options)
            D = min(D, D1)
            [tmp, i] = max(D(: ))
            [x, y] = ind2sub([n n], i)
            vertex(: , end + 1) = [x; y]
        % compute the Delaunay triangulation
        [D, Z, Q] = perform_fast_marching(1./ W, vertex)
        faces = compute_voronoi_triangulation(Q, vertex)
        % compute approximation
        vgeod = compute_orthoproj_triangulation(vertex, faces, M)
        Mgeod = compute_triangulation_interpolation(faces, vertex, vgeod, n)
        % display
        imageplot(clamp(Mgeod), ['\alpha = ' num2str(alpha) ', SNR = ' num2str(snr(Mgeod, M), 3) 'dB'], 2, 2, ialpha)


