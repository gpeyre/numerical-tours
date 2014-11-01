def exo1():
    """
    Compute curves joining the start point to several points along the
    boundary.
    """
    ms = 20
    U = geod(x)
    npaths = 30
    sel = round(linspace(1, L + 1, npaths + 1)); sel(end) = []
    end_points = b(: , sel)
    clf; hold on
    imageplot(1-S)
    for i in 1: npaths:
        p = compute_geodesic(U, end_points(: , i))
        h = plot(p(2, : ), p(1, : ), 'g'); set(h, 'LineWidth', lw)
    h = plot(x(2), x(1), '.r'); set(h, 'MarkerSize', ms)
    h = plot(end_points(2, : ), end_points(1, : ), '.b'); set(h, 'MarkerSize', ms)
    axis ij


def exo2():
    """
    Compute the geodesic distance matrix $\de$.
    """
    delta = zeros(N, N)
    sel = X(1, : ) + (X(2, : )-1)*q
    for i in 1: N:
        U = geod(X(: , i))
        delta(: , i) = U(sel)
    delta = (delta + delta')/ 2


def exo3():
    """
    Perform the SMACOF iterative algorithm.
    Save in a variable |s(l)| the values of
    Stress$( X^{(\ell)} )$.
    """
    niter = 50
    s = []
    Y = X/ q
    ndisp = [1 5 niter Inf]
    clf; k = 1
    hold on
    Y = Y-repmat(mean(Y, 2), [1 N])
    h = plot(Y(2, [1: N0 1]), Y(1, [1: N0 1]))
    for i in 1: niter:
        Y = Y * B(D(Y))' / N
        % update
        Y = Y-repmat(mean(Y, 2), [1 N])
        % record stress
        s(end + 1) = Stress(D(Y))
        if ndisp(k) = =i
            plot(Y(2, [1: N0 1]), Y(1, [1: N0 1]))
            axis('equal'); axis('off')
            k = k + 1
    axis('equal'); axis('off'); axis('ij')


def exo4():
    """
    Implement a shape retrival algorithm based on these bending invariants.
    o correction for this exercise.
    """


