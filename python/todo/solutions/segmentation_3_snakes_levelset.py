def exo1():
    """
    Load a square shape $\phi_2$ at a different position for the center.
    adius
    enter
    hape
    """
    do_fast = 1
    r = n/ 3
    c = n  - 10 - [r r]
    phi2 = max(abs(X-c(1)), abs(Y-c(2))) - r


def exo2():
    """
    Compute the intersection and the union of the two shapes.
    Store the union in $\phi_0$ (variable |phi0|),
    that we will use in the remaining part of the tour.
    """
    phi0 = min(phi1, phi2)
    subplot(1, 2, 1)
    plot_levelset(min(phi1, phi2))
    title('Union')
    subplot(1, 2, 2)
    plot_levelset(max(phi1, phi2))
    title('Intersection')


def exo3():
    """
    Implement the mean curvature motion.
    """
    if do_fast
    phi = phi0; % initialization
    k = 0
    for i in 1: niter:
        g0 = grad(phi, options)
        d = max(eps, sqrt(sum(g0.^2, 3)))
        g = g0 ./ repmat(d, [1 1 2])
        K = d .* div(g, options)
        phi = phi + tau*K
        % redistance the function from time to time
        if mod(i, 30) = =0
            % phi = perform_redistancing(phi)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            subplot(2, 2, k)
            plot_levelset(phi)


def exo4():
    """
    Compute an initial shape $\phi_0$ at time $t=0$,
    for instance a centered square.
    """
    [Y, X] = meshgrid(1: n, 1: n)
    r = n/ 3
    c = [n n]/ 2
    phi0 = max(abs(X-c(1)), abs(Y-c(2))) - r


def exo5():
    """
    Compute and store in |G| the gradient $G(\phi)$ (right hand side of the PDE)
    using the current value of the distance function $\phi$.
    ormalized gradient
    radient
    """
    gD = grad(phi, options)
    d = max(eps, sqrt(sum(gD.^2, 3)))
    g = gD ./ repmat(d, [1 1 2])
    G = - W .* d .* div(g, options) - sum(gW.*gD, 3)


def exo6():
    """
    Implement the geodesic active contours gradient descent.
    Do not forget to do the re-distancing.
    """
    if do_fast
    phi = phi0
    k = 0
    gW = grad(W, options)
    for i in 1: niter:
        gD = grad(phi, options)
        d = max(eps, sqrt(sum(gD.^2, 3)))
        g = gD ./ repmat(d, [1 1 2])
        G = W .* d .* div(g, options) + sum(gW.*gD, 3)
        phi = phi + tau*G
        if mod(i, 30) = =0
            phi = perform_redistancing(phi)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            subplot(2, 2, k)
            plot_levelset(phi, 0, f0)


def exo7():
    """
    Compute an initial level set function $\phi_0$, stored in |phi0|,
    for instance many small circles.
    """
    [Y, X] = meshgrid(1: n, 1: n)
    k = 4; %number of circles
    r = .3*n/ k
    phi0 = zeros(n, n) + Inf
    for i in 1: k:
    for j in 1: k:
            c = ([i j]-1)*(n/ k) + (n/ k)*.5
            phi0 = min(phi0, sqrt((X-c(1)).^2 + (Y-c(2)).^2) - r)
    subplot(1, 2, 1)
    plot_levelset(phi0)
    subplot(1, 2, 2)
    plot_levelset(phi0, 0, f0)


def exo8():
    """
    Compute this gradient $G(\phi)$ using the current value of the distance function
    (phi$.
    radient
    """
    gD = grad(phi, options)
    d = max(eps, sqrt(sum(gD.^2, 3)))
    g = gD ./ repmat(d, [1 1 2])
    G = d .* div(g, options) - lambda*(f0-c1).^2 + lambda*(f0-c2).^2


def exo9():
    """
    Implement the full gradient descent.
    """
    if do_fast
    phi = phi0
    k = 0
    for i in 1: niter:
        gD = grad(phi, options)
        d = max(eps, sqrt(sum(gD.^2, 3)))
        g = gD ./ repmat(d, [1 1 2])
        G = d .* div(g, options) - lambda*(f0-c1).^2 + lambda*(f0-c2).^2
        phi = phi + tau*G
        if mod(i, 30) = =0
            phi = perform_redistancing(phi)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            subplot(2, 2, k)
            plot_levelset(phi, 0, f0)


def exo10():
    """
    In the case that one does not know precisely the constants $c_1$ and $c_2$,
    how to update them automatically during the evolution ? Implement this method.
    """


