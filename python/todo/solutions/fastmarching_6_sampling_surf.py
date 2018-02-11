def exo1():
    """
    Perform the farthest point sampling of |m=500| points.
    nitialize
    """
    m = 500
    landmarks = [100]
    [D, Z, Q] = perform_fast_marching_mesh(vertex, faces, landmarks)
    k = 1; displist = [10 50 100 m]
    for i in 2: m:
        % select
        [tmp, landmarks(end + 1)] = max(D)
        % update
        options.constraint_map = D
        [D1, Z, Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options)
        D = min(D, D1)
        if i = =displist(k)
            subplot(2, 2, k)
            hold on
            options.face_vertex_color = perform_hist_eq(D, 'linear')
            plot_mesh(vertex, faces, options)
            colormap jet(256)
            h = plot3(vertex(1, landmarks), vertex(2, landmarks), vertex(3, landmarks), 'r.')
            set(h, 'MarkerSize', 20)
            k = k + 1


def exo2():
    """
    Perform a spacially adative remeshing.
    nitialize
    """
    m = 200
    landmarks = [100]
    [D, Z, Q] = perform_fast_marching_mesh(vertex, faces, landmarks)
    k = 1; displist = [100 m]
    for i in 2: m:
        % select
        [tmp, landmarks(end + 1)] = max(D)
        % update
        options.constraint_map = D
        [D1, Z, Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options)
        D = min(D, D1)
        if i = =displist(k)
            [D, Z, Q] = perform_fast_marching_mesh(vertex, faces, landmarks)
            % compute the mesh
            V = Q(faces); V = sort(V, 1)
            V = unique(V', 'rows')'
            d = 1 + (V(1, : )~ = V(2, : )) + (V(2, : )~ = V(3, : ))
            %
            I = find(d = =3); I = sort(I)
            z = zeros(n, 1)
            z(landmarks) = (1: length(landmarks))'
            facesV = z(V(: , I))
            vertexV = vertex(: , landmarks)
            % Re-orient the faces so that they point outward of the mesh.
            options.method = 'slow'
            options.verb = 0
            facesV = perform_faces_reorientation(vertexV, facesV, options)
            % display
            subplot(1, 2, k)
            options.face_vertex_color = []
            plot_mesh(vertexV, facesV, options)
            shading faceted
            %
            k = k + 1


def exo3():
    """
    Design a metric |W| so that the sampling is densed in area where |C| is
    large.
    
    isplay
    """
    W = rescale(min(C, .1), .001, 1)
    options.W = W
    landmarks = [5000]
    options.constraint_map = []
    [D, Z, Q] = perform_fast_marching_mesh(vertex, faces, landmarks, options)
    hold on
    options.face_vertex_color = mod(20*D/ max(D), 1)
    plot_mesh(vertex, faces, options)
    colormap jet(256)
    h = plot3(vertex(1, landmarks), vertex(2, landmarks), vertex(3, landmarks), 'r.')
    set(h, 'MarkerSize', 20)


def exo4():
    """
    Use such a metric to perform feature sensitive remeshing.
    Tune the metric to reduce as much as possible the Hausdorff
    approximation error.
    """


