def exo1():
    """
    Using |Q|, compute the faces |faces| of the Delaunay triangulation.
    To that end,
    extract each quad of values |Q(i,j),Q(i+1,j),Q(i+1,j+1),Q(i,j+1)|,
    and add a new face when three of these four values are different (this
    corresponds to a Voronoi point).
    Display the obtained triangulation.
    ompute quad of values
    """
    V = []
    v = Q(1: end-1, 1: end-1); V = [V v(: )]
    v = Q(2: end, 1: end-1); V = [V v(: )]
    v = Q(1: end-1, 2: end); V = [V v(: )]
    v = Q(2: end, 2: end); V = [V v(: )]
    V = sort(V, 2)
    V = unique(V, 'rows')
    d = (V(: , 1)~ = V(: , 2)) + (V(: , 2)~ = V(: , 3)) + (V(: , 3)~ = V(: , 4))
    V = V(d = =2, : )
    for i in 1: size(V, 1):
        V(i, 1: 3) = unique(V(i, : ))
    faces = V(: , 1: 3)'


def exo2():
    """
    Iterate the sampling process to add more and more points.
    """
    nmax = 200
    ndisp = [20 50 100 200]
    vertex = [1; 1]
    [D, Z, Q] = perform_fast_marching(1./ W, vertex)
    clf; u = 1
    vertex_svg = {}
    faces_svg = {}
    for k in 2: nmax:
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
            u = u + 1
            % compute the Delaunay triangulation
            [D1, Z, Q] = perform_fast_marching(1./ W, vertex)
            vertex_svg{end + 1} = vertex
            faces_svg{end + 1} = compute_voronoi_triangulation(Q, vertex)


def exo3():
    """
    Display the geodesic Delaunay triangulation corresponding to the
    sampling
    """
    for i in 1: 4:
        subplot(2, 2, i)
        plot_triangulation(vertex_svg{i}, faces_svg{i}, M)
    colormap jet(256)


def exo4():
    """
    Perform the Lloyd iterative algorithm.
    etric
    eed random points.
    """
    n = 512; W = ones(n)
    p = 40
    vertex = floor(rand(2, p)*(n-1)) + 1
    disp_list = [1, 2, 5, 20]; k = 1
    clf; hold on
    for i in 1: max(disp_list):
        % Compute Voronoi partition.
        [D, Z, Q] = perform_fast_marching(1./ W, vertex)
        if i = =disp_list(k)
            subplot(2, 2, k)
            imageplot(Q'); hold on
            h = plot(vertex(1, : ), vertex(2, : ), 'k.')
            set(h, 'MarkerSize', 15)
            colormap(jet(256))
            k = k + 1
        % Re-center each point at the barycenter of its cell.    
    for i in 1: p:
            [x, y] = ind2sub(size(W), find(Q = =i))
            vertex(: , i) = [mean(x); mean(y)]


def exo5():
    """
    Perform the Lloyd iterative algorithm.
    eed random points.
    """
    p = 150
    vertex = floor(rand(2, p)*(n-1)) + 1
    disp_list = [1, 2, 5, 30]; k = 1
    clf; hold on
    for i in 1: max(disp_list):
        % Compute Voronoi partition.
        [D, Z, Q] = perform_fast_marching(1./ W, vertex)
        if i = =disp_list(k)
            subplot(2, 2, k)
            imageplot(Q'); hold on
            h = plot(vertex(1, : ), vertex(2, : ), 'k.')
            set(h, 'MarkerSize', 15)
            colormap(jet(256))
            k = k + 1
        % Re-center each point at the barycenter of its cell.    
    for i in 1: p:
            [x, y] = ind2sub(size(W), find(Q = =i))
            vertex(: , i) = [mean(x); mean(y)]


