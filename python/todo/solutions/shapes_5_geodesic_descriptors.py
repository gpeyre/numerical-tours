def exo1():
    """
    Using |options.nb_iter_max| display the progression of the Fast
    Marching.
    """
    p = sum(M(: ) = =1)
    qlist = ceil([.1 .3 .7 1]*p)
    for i in 1: 4:
        options.nb_iter_max = qlist(i)
        D = perform_fast_marching(W, start_points, options)
        subplot(2, 2, i)
        display_shape_function(perform_hist_eq(D, 'linear'))
        hold on
        h = plot(bound(2, : ), bound(1, : ), 'k')
        set(h, 'LineWidth', 2)
        axis('ij')
    options.nb_iter_max = Inf


def exo2():
    """
    Compute curves joining the start point to several points along the
    boundary.
    """
    npaths = 30
    sel = round(linspace(1, nbound + 1, npaths + 1)); sel(end) = []
    end_points = bound(: , sel)
    clf; hold on
    imageplot(1-M)
    for i in 1: npaths:
        p = compute_geodesic(D, end_points(: , i))
        h = plot(p(2, : ), p(1, : ), 'g'); set(h, 'LineWidth', lw)
    h = plot(start_points(2), start_points(1), '.r'); set(h, 'MarkerSize', ms)
    h = plot(end_points(2, : ), end_points(1, : ), '.b'); set(h, 'MarkerSize', ms)
    axis ij


def exo3():
    """
    Build a collection |E| of distance maps, so that |E(:,:,i)| is the
    geodesic distance to |samples(:,i)|.
    """
    E = zeros(n, n, nb_samples)
    for i in 1: nb_samples:
        % progressbar(i, nb_samples)
        [d, tmp] = perform_fast_marching(W, samples(: , i), options)
        d(M = =0) = 0
        d(d = =Inf) = 0
        d(d >1e5) = 0
        E(: , : , i) = d


def exo4():
    """
    Load a library of shapes. Compute the different histograms for these
    shapes.
    """


def exo5():
    """
    Perform the retrieval by comparing the histogram. Test diffetent metrics
    for the retrieval.
    """


