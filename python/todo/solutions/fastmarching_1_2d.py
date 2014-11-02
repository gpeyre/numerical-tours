def exo1():
    """
    Using |options.nb_iter_max|, display the progressive propagation.
    This corresponds to displaying the front
    $ \enscond{x}{D(x) \leq t} $ for various arrival times $t$.
    """
    niter = round(linspace(.1, 1, 6)*n^2)
    for i in 1: length(niter):
        options.nb_iter_max = niter(i)
        options.end_points = []
        [D, S] = perform_fast_marching(1./ W, x0, options)
        subplot(2, 3, i)
        hold on
        imageplot(convert_distance_color(D, f))
        h = plot(x0(2), x0(1), '.r'); set(h, 'MarkerSize', 25)


def exo2():
    """
    Perform the full geodesic path extraction by iterating the gradient
    descent. You must be very careful when the path become close to
    $x_0$, because the distance function is not differentiable at this
    point. You must stop the iteration when the path is close to $x_0$.
    """
    gamma = x1
    for i in 1: 1.5*n/ tau:
        gamma(: , end + 1) = gamma(: , end) - tau*Geval(G, gamma(: , end))
        if norm(gamma(: , end)-x0) <1
            break
    gamma(: , end + 1) = x0


def exo3():
    """
    Study the influence of the $\epsilon$ parameter.
    """
    elist = [.1 .5 1 10]
    for ie in 1: length(elist):
        epsilon = elist(ie)
        W = epsilon + abs(f-f(x0(1), x0(2)))
        [D, S] = perform_fast_marching(1./ W, x0, options)
        G0 = grad(D, options)
        G = G0 ./ repmat(sqrt(sum(G0.^2, 3)), [1 1 2])
        % 
        gamma = x1
    for i in 1: 1.5*n/ tau:
            gamma(: , end + 1) = gamma(: , end) - tau*Geval(G, gamma(: , end))
            if norm(gamma(: , end)-x0) <1
                break
        gamma(: , end + 1) = x0
        %
        subplot(2, 2, ie); hold on
        imageplot(f)
        h = plot(gamma(2, : ), gamma(1, : ), '.b'); set(h, 'LineWidth', 2)
        h = plot(x0(2), x0(1), '.r'); set(h, 'MarkerSize', 25)
        h = plot(x1(2), x1(1), '.b'); set(h, 'MarkerSize', 25)
        axis ij
        title(['\epsilon = ' num2str(epsilon)])


def exo4():
    """
    Perform the shortest path
    extraction for various images such as 'cavern' or 'mountain'.
    oad
    radient
    isplay
    """
    f = load_image('cavern', n)
    epsilon = 1e-2
    W = epsilon + rescale(-f)
    x0 = [45; 280]; x1 = [275; 25]
    options.nb_iter_max = Inf
    options.end_points = []; % x1
    [D, S] = perform_fast_marching(1./ W, x0, options)
    G0 = grad(D, options)
    G = G0 ./ repmat(sqrt(sum(G0.^2, 3)), [1 1 2])
    gamma = x1
    for i in 1: 1.5*n/ tau:
        gamma(: , end + 1) = gamma(: , end) - tau*Geval(G, gamma(: , end))
        if norm(gamma(: , end)-x0) <1
            break
    gamma(: , end + 1) = x0
    subplot(1, 2, 1)
    hold on
    imageplot(repmat(f, [1 1 3]), 'Image')
    h = plot(gamma(2, : ), gamma(1, : ), '.b'); set(h, 'LineWidth', 2)
    h = plot(x0(2), x0(1), '.r'); set(h, 'MarkerSize', 25)
    h = plot(x1(2), x1(1), '.b'); set(h, 'MarkerSize', 25)
    subplot(1, 2, 2)
    hold on
    imageplot(convert_distance_color(perform_hist_eq(D, 'linear'), f), 'Distance')
    h = plot(gamma(2, : ), gamma(1, : ), '.k'); set(h, 'LineWidth', 2)
    h = plot(x0(2), x0(1), '.r'); set(h, 'MarkerSize', 25)
    h = plot(x1(2), x1(1), '.b'); set(h, 'MarkerSize', 25)


def exo5():
    """
    Extract the set of points that are along the boundary of the Voronoi
    region. This corresponds for instance to the points of the region
    $ \enscond{x}{Q(x)=1} $
    that have one neighbor inside the region
    $ \enscond{x}{Q(x)=2} $.
    Compute the geodesic distance $D(x)$ at these points, and choose two points
    $a$ and $b$ on this boundary that have small values of $D$.
    int: you can use a convolution |U=conv2(double(Q==2),h,'same')| with a
    ell chose kernel |h| to located the points |U>0| with at least 1
    eighbor.
    
    ubplot(2,1,1);
    
    """
    h = [0 1 0; 1 0 1; 0 1 0]
    B = (Q = =1) & (conv2(double(Q = =2), h, 'same') >0)
    U = find(B)
    [xa, xb] = ind2sub(size(f), U)
    [xa, I] = sort(xa); xb = xb(I); U = U(I)
    dU = D(U)
    k = [65 259]
    x1 = [[xa(k(1)); xb(k(1))] [xa(k(2)); xb(k(2))]]
    hold on
    imageplot(A, 'Boundary points')
    h = plot(x0(2, : ), x0(1, : ), '.g'); set(h, 'MarkerSize', 25)
    h = plot(xb, xa, 'g'); set(h, 'LineWidth', 2)
    h = plot(x1(2, : ), x1(1, : ), '.b'); set(h, 'MarkerSize', 25)
    if 0
    subplot(2, 1, 2)
    hold on
    h = plot([k(1) k(1)], [min(dU) max(dU)], 'r: '); set(h, 'LineWidth', 2)
    h = plot([k(2) k(2)], [min(dU) max(dU)], 'r: '); set(h, 'LineWidth', 2)
    h = plot(dU); axis('tight'); title('D along the boundary'); set(h, 'LineWidth', 2)


def exo6():
    """
    Extract the geodesics joining $a$ and $b$ to the two starting points
    (this makes 4 geodesic curves). Use them to perform segmentation.
     D1 = D; D1(D1==Inf) = max(D1(D1~=Inf));
    isplay the curves
    """
    options.nb_iter_max = Inf
    options.end_points = []
    tau = .8
    curves = {}
    for i in 1: 2:
        % FM
        [D, S, Q] = perform_fast_marching(1./ W, x0(: , i), options)
        % gradient
        G = grad(D, options)
        G = G ./ repmat(sqrt(sum(G.^2, 3) + 1e-9), [1 1 2])
        % Geodesic
    for j in 1: 2:
            gamma = x1(: , j)
    for iter in 1: 1.5*n/ tau:
                g = [interp2(1: n, 1: n, G(: , : , 1), gamma(2, end), gamma(1, end)); ...
                    interp2(1: n, 1: n, G(: , : , 2), gamma(2, end), gamma(1, end))]
                gamma(: , end + 1) = gamma(: , end) - tau*g
                if norm(gamma(: , end)-x0(: , i)) <1
                    break
            gamma(: , end + 1) = x0(: , i)
            curves{end + 1} = gamma
    col = {'r', 'g', 'b', 'c'}
    hold on
    imageplot(A, 'Boundary points')
    for i in 1: length(curves):
        c = curves{i}
        h = plot(c(2, : ), c(1, : ), col{i}); set(h, 'LineWidth', 2)
    h = plot(x0(2, : ), x0(1, : ), '.g'); set(h, 'MarkerSize', 25)
    h = plot(x1(2, : ), x1(1, : ), '.b'); set(h, 'MarkerSize', 25)


def exo7():
    """
    Perform partial propagations from $x_0$.
    """
    niter = round(linspace(.1, 1, 6)*n^2)
    for i in 1: length(niter):
        options.nb_iter_max = niter(i)
        options.end_points = []
        [D, S] = perform_fast_marching(1./ W, x0, options)
        subplot(2, 3, i)
        hold on
        imageplot(convert_distance_color(D, f))
        h = plot(x0(2), x0(1), '.r'); set(h, 'MarkerSize', 25)


def exo8():
    """
    Extract geodesics joining several points $x_1$ to the central point
    $x_0$.
    radient
    xtract centerlines
    isplay the curves
    """
    x1 = [[175; 5] [21; 42] [48; 133] [244; 78] [191; 40] ...
             [100; 13] [66; 42] [183; 66] [220; 117]]
    D1 = D; D1(D1 = =Inf) = max(D1(D1~ = Inf))
    G = grad(D1, options)
    G = G ./ repmat(sqrt(sum(G.^2, 3) + 1e-9), [1 1 2])
    curves = {}
    for k in 1: size(x1, 2):
        % extract curve
        gamma = x1(: , k)
    for iter in 1: 1.5*n/ tau:
            g = [interp2(1: n, 1: n, G(: , : , 1), gamma(2, end), gamma(1, end)); ...
                interp2(1: n, 1: n, G(: , : , 2), gamma(2, end), gamma(1, end))]
            gamma(: , end + 1) = clamp(gamma(: , end) - tau*g, 1, n)
            if norm(gamma(: , end)-x0) <1
                break
        gamma(: , end + 1) = x0
        curves{end + 1} = gamma
    hold on
    imageplot(f, 'Boundary points')
    for i in 1: length(curves):
        c = curves{i}
        h = plot(c(2, : ), c(1, : ), 'r'); set(h, 'LineWidth', 2)
    h = plot(x0(2, : ), x0(1, : ), '.g'); set(h, 'MarkerSize', 25)
    h = plot(x1(2, : ), x1(1, : ), '.b'); set(h, 'MarkerSize', 25)


def exo9():
    """
    Perform the dual propagation, and stop it when the front meet.
    Extract the two half geodesic curves.
    ual propagation.
    xtract first the geodesic paths
    terations
    """
    options.end_points = []
    iterlist = .37*[.25 .5 .75 1]*n^2
    options.nb_iter_max = Inf
    [D, S] = perform_fast_marching(1./ W, x0(: , 1), options)
    gamma = compute_geodesic(D, x0(: , 2))
    for i in 1: 4:
        options.nb_iter_max = iterlist(i)
        [D, S] = perform_fast_marching(1./ W, x0, options)
        subplot(2, 2, i)
        hold on
        imageplot(convert_distance_color(D, f))
        if i = =4
            h = plot(gamma(2, : ), gamma(1, : ), 'k'); set(h, 'LineWidth', 2)
            % select extremal point
            u = interp2(1: n, 1: n, D, gamma(2, : ), gamma(1, : ))
            [tmp, i] = max(u); q = gamma(: , i)
            h = plot(q(2, : ), q(1, : ), '.b'); set(h, 'MarkerSize', 25)
        h = plot(x0(2, : ), x0(1, : ), '.r'); set(h, 'MarkerSize', 25)


