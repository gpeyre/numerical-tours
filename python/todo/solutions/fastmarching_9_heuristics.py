def exo1():
    """
    Display the set of points satisfying |D+H<=T| for several value of the
    threshold |T>=D(pend)|. What do you observe ?
    """
    d = (D(pend(1), pend(2)) + H(pstart(1), pstart(2)))/ 2
    Tlist = [1.01 1.05 1.1 1.2] * d
    for i in 1: 4:
        t = Tlist(i)
        U = cat(3, M, M, M)
        I = find(D + H <= t)
        U(I) = 1; U([I + n^2, I + 2*n^2]) = U([I + n^2, I + 2*n^2])*.3
        subplot(2, 2, i)
        hold on
        imageplot(U)
        h = plot(p(2, : ), p(1, : ), '.k'); set(h, 'LineWidth', 2)
        h = plot(pstart(2), pstart(1), '.g'); set(h, 'MarkerSize', 25)
        h = plot(pend(2), pend(1), '.b'); set(h, 'MarkerSize', 25)
        axis ij


def exo2():
    """
    Display the explored region for different values of |weight|.
    """
    wlist = [.6 .8 .9 .95]
    for i in 1: 4:
        weight = wlist(i)
        options.end_points = pend
        options.heuristic = weight*H
        options.nb_iter_max = Inf
        options.constraint_map = Inf + zeros(n)
        [D, S] = perform_fast_marching(1./ W, pstart, options)
        %
        I = find(S <0)
        U = cat(3, M, M, M)
        U(I) = 1; U([I + n^2, I + 2*n^2]) = U([I + n^2, I + 2*n^2])*.3
        subplot(2, 2, i)
        hold on
        imageplot(U)
        h = plot(p(2, : ), p(1, : ), '.k'); set(h, 'LineWidth', 2)
        h = plot(pstart(2), pstart(1), '.g'); set(h, 'MarkerSize', 25)
        h = plot(pend(2), pend(1), '.b'); set(h, 'MarkerSize', 25)
        axis ij


def exo3():
    """
    Display the convergence of the heuristic as the number of landmark
    increases.
    """
    qlist = [2 5 10 50]
    q = max(qlist)
    landmarks = floor(rand(2, q)*n) + 1
    Dland = zeros(n, n, q)
    for i in 1: q:
        Dland(: , : , i) = perform_fast_marching(1./ W, landmarks(: , i))
    for i in 1: 4:
        q = qlist(i)
        Dend = Dland(pend(1), pend(2), : )
        H = max(abs(Dland(: , : , 1: q)-repmat(Dend(1: q), [n n 1])), [], 3)
        subplot(2, 2, i)
        hold on
        imageplot(H)
        contour(H, 10, 'k', 'LineWidth', 2)
        colormap jet(256)
        h = plot(landmarks(1, 1: q), landmarks(2, 1: q), 'y.')
        set(h, 'MarkerSize', 15)
        axis ij


def exo4():
    """
    Perform the heuristically driven propagation with a landmark-based
    heuristic.
    """
    for i in 1: 4:
        q = qlist(i)
        Dend = Dland(pend(1), pend(2), : )
        H = max(abs(Dland(: , : , 1: q)-repmat(Dend(1: q), [n n 1])), [], 3)
        %
        options.end_points = pend
        options.heuristic = H
        options.nb_iter_max = Inf
        options.constraint_map = Inf + zeros(n)
        [D, S] = perform_fast_marching(1./ W, pstart, options)
        %
        I = find(S <0)
        U = cat(3, M, M, M)
        U(I) = 1; U([I + n^2, I + 2*n^2]) = U([I + n^2, I + 2*n^2])*.3
        subplot(2, 2, i)
        hold on
        imageplot(U)
        h = plot(p(2, : ), p(1, : ), '.k'); set(h, 'LineWidth', 2)
        h = plot(pstart(2), pstart(1), '.g'); set(h, 'MarkerSize', 25)
        h = plot(pend(2), pend(1), '.b'); set(h, 'MarkerSize', 25)
        h = plot(landmarks(1, 1: q), landmarks(2, 1: q), 'y.'); set(h, 'MarkerSize', 15)
        axis ij


def exo5():
    """
    Find a strategy to find optimal seeding position for the landmarks.
    """


def exo6():
    """
    Perform the landmark-based heuristically driven propagation on a 3D
    mesh.
    """


