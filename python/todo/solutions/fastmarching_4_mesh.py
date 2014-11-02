def exo1():
    """
    Using |options.nb_iter_max|, display the progression of the propagation.
    """
    nblist = round(linspace(.1, 1, 6)*nvert)
    for i  in  1: length(nblist):
        options.nb_iter_max = nblist(i)
        [D, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
        subplot(2, 3, i)
        col = D; col(col = =Inf) = 0
        options.face_vertex_color = col
        hold('on')
        plot_mesh(vertex, faces, options)
        colormap jet(256)
        % display here starting points


def exo2():
    """
    For each point |pend(k)|, compute a discrete geodesic path |path| such
    that |path(1)=pend(k)| and |D(path(i+1))<D(path(i))|
    with |[path(i), path(i+1)]| being an edge of the mesh.
    This means that |path(i+1)| is an element of |vring{path(i)}|.
    Display the paths on the mesh.
    """
    pathsD = {}
    for k in 1: nend:
        % path purely on edges
        vlist = pend(k)
        vprev = D(vlist(end))
        while true
            x0 = vlist(end)
            r = vring{x0}
            [v, J] = min(D(r))
            x = r(J)
            if v >= vprev || v = =0
                break
            vprev = v
            vlist(end + 1) = x
        pathsD{end + 1} = vertex(: , vlist)
    plot_fast_marching_mesh(vertex, faces, Q, pathsD, options)


def exo3():
    """
    Using |options.nb_iter_max|, display the progression of the propagation for constant |W|.
    """
    options.W = ones(nvert, 1)
    nblist = round([.05 .15 .4 .6]*nvert)
    for i  in  1: length(nblist):
        options.nb_iter_max = nblist(i)
        [D0, S, Q0] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
        subplot(2, 2, i)
        plot_fast_marching_mesh(vertex, faces, D0, [], options)


def exo4():
    """
    Using |options.nb_iter_max|, display the progression of the propagation for a curvature based |W|.
    """
    options.W = W
    for i  in  1: length(nblist):
        options.nb_iter_max = nblist(i)
        [D, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
        subplot(2, 2, i)
        plot_fast_marching_mesh(vertex, faces, D, [], options)


def exo5():
    """
    Extract geodesics.
    ompute distances
    ompute paths
    isplay
    """
    options.nb_iter_max = Inf
    options.W = ones(nvert, 1)
    [D0, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
    options.W = W
    [D, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
    pend = [7175; 17991; 2293]
    options.method = 'continuous'
    paths0 = compute_geodesic_mesh(D0, vertex, faces, pend, options)
    paths = compute_geodesic_mesh(D, vertex, faces, pend, options)
    subplot(1, 2, 1)
    plot_fast_marching_mesh(vertex, faces, perform_hist_eq(D0, 'linear'), paths0, options)
    subplot(1, 2, 2)
    plot_fast_marching_mesh(vertex, faces, perform_hist_eq(D, 'linear'), paths, options)


def exo6():
    """
    Using |options.nb_iter_max|, display the progression of the propagation for a curvature based |W|.
    """
    options.W = W
    nblist = round([.05 .15 .4 .6]*nvert)
    for i  in  1: length(nblist):
        options.nb_iter_max = nblist(i)
        [D, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
        subplot(2, 2, i)
        plot_fast_marching_mesh(vertex, faces, D, [], options)


def exo7():
    """
    Extract geodesics.
    ompute distances
    ompute paths
    isplay
    """
    options.nb_iter_max = Inf
    options.W = ones(nvert, 1)
    [D0, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
    options.W = W
    [D, S, Q] = perform_fast_marching_mesh(vertex, faces, pstarts, options)
    pend = [9969; 5073]
    options.method = 'continuous'
    paths0 = compute_geodesic_mesh(D0, vertex, faces, pend, options)
    paths = compute_geodesic_mesh(D, vertex, faces, pend, options)
    subplot(1, 2, 1)
    plot_fast_marching_mesh(vertex, faces, f, paths0, options)
    subplot(1, 2, 2)
    plot_fast_marching_mesh(vertex, faces, f, paths, options)
    colormap gray(256)


