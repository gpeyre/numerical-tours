def exo1():
    """
    Perform the full iterative algorithm until convergence.
    Store in |err(l)| the fixed point error
    $ \norm{ U^{(\ell+1)} - U^{(\ell)} } $.
    _Note:_ you might want to put outside of the loop
    all the quantities that do not depend on $u$, e.g.
    $S, a, $ etc.
    isplay
    """
    options.niter = 250*3
    options.verb = 0
    options.svg_rate = 10
    options.U = []
    [U, err, Usvg] = perform_geodesic_iterative(vertex, faces, W, I, options)
    t = round(linspace(1, size(Usvg, 2)*.5, 4))
    for kdisp in 1: 4:
        subplot(2, 2, kdisp)
        options.face_vertex_color = Usvg(: , t(kdisp))
        plot_mesh(vertex, faces, options)
        colormap jet(256)
        shading interp


def exo2():
    """
    Compute and display the geodesic distance.
    """
    options.niter = 200
    options.U = sqrt(sum(vertex.^2))'
    [U, err] = perform_geodesic_iterative(vertex, faces, W, I, options)
    options.face_vertex_color = mycolor(U, 3)
    plot_mesh(vertex, faces, options)
    colormap jet(256); shading interp


def exo3():
    """
    Display the convergence of the computed geodesic distance to the the true
    geodesic distance (which is the Euclidean distance $ \norm{x_i} $) as $n$
    increases.
    _Note:_ the triangulation with increasing number of points should be
    refining (i.e. a finer triangulation should contains all the other
    ones).
    """
    nlist = round(linspace(100, 4000, 10))
    options.niter = 200
    err = []
    for i in 1: length(nlist):
        n = nlist(i)
        warning off; rand('state', 0); warning on
        vertex = [2*rand(2, n)-1; zeros(1, n)]
        vertex(: , 1) = 0
        faces = compute_delaunay(vertex)
        options.U = sqrt(sum(vertex.^2))'
        [U, e] = perform_geodesic_iterative(vertex, faces, W, I, options)
        err(end + 1) = norm(U-options.U)/ sqrt(n)
    h = plot(nlist, err, '.-')
    set(h, 'LineWidth', 2)
    axis tight


def exo4():
    """
    Compute the geodesic distance for a metric $W_i$ that is not constant
    over the mesh.
    
    """
    q = 12
    vertex(1, 2: q + 1) = .97
    vertex(2, 2: q + 1) = linspace(0.03, .97, q)
    faces = compute_delaunay(vertex)
    W = ones(n, 1);  W(vertex(1, : ) <.5) = 1/ 2
    I = 1
    options.niter = 300
    options.U = sqrt(sum((vertex-repmat(vertex(: , I), [1 n])).^2))'
    [U, err] = perform_geodesic_iterative(vertex, faces, W, I, options)
    options.method = 'continuous'
    paths = compute_geodesic_mesh(U, vertex, faces, 2: q + 1, options)
    plot_fast_marching_mesh(vertex, faces, mycolor(U, 8), paths, options)


