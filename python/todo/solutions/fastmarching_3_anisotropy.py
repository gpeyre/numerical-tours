def exo1():
    """
    Compute the geodesic distance for several anisotropy, and for several
    starting points.
    """
    npoints = 10
    nstart = floor(rand(2, npoints)*n) + 1
    aniso_list = [.02 .1 .2 .5]
    for i in 1: length(aniso_list):
        anisotropy = aniso_list(i)
        H = perform_tensor_recomp(e1, e2, ones(n), ones(n)*1/ anisotropy)
        [D, dUx, dUy, Vor, L] = fm2dAniso([hx; hy], H, pstart)
        imageplot(convert_distance_color(D), ['Anisotropy = ' num2str(anisotropy)], 2, 2, i)


def exo2():
    """
    Perform farthest point sampling.
    """
    nmax = 200
    ndisp = [20 50 100 200]
    vertex = [1; 1]
    [D, dUx, dUy, Vor, L] = fm2dAniso([hx; hy], H, vertex)
    clf; u = 1
    vertex_svg = {}
    faces_svg = {}
    for k in 2: nmax:
        options.constraint_map = D
        [D, dUx, dUy, Vor, L] = fm2dAniso([hx; hy], H, vertex)
        [tmp, i] = max(D(: ))
        [x, y] = ind2sub([n n], i)
        vertex(: , end + 1) = [x; y]
        if k = =ndisp(u)
            subplot(2, 2, u)
            hold on
            imageplot(D, [num2str(k) ' points'])
            plot(vertex(2, : ), vertex(1, : ), 'r.')
            u = u + 1
            % compute the Delaunay triangulation
            % [D1, Z, Q] = perform_fast_marching(1./ W, vertex)
            % vertex_svg{end + 1} = vertex
            % faces_svg{end + 1} = compute_voronoi_triangulation(Q, vertex)


def exo3():
    """
    Compute a metric |H| adapted to the approximation of this image.
    """


def exo4():
    """
    Perform farthest point sampling.
    """


def exo5():
    """
    Compute the geodesic Delaunay triangulation of this set of point.
    """


def exo6():
    """
    Perform image approximation using linear splines.
    """


