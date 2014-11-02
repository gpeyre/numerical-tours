def exo1():
    """
    Compute a set of |n = 300| vertex by iterating this farthest
    point sampling. Display the progression of the sampling.
    """
    n = 400
    landmarks = 1
    Dland = []
    k = 1
    displ = round(linspace(0, 1, 5)*n); displ(1) = []
    for i in 1: n:
        if not(isempty(Dland))
            [tmp, landmarks(end + 1)] = max(min(Dland, [], 2))
        [Dland(: , end + 1), S, Q] = perform_fast_marching_mesh(vertex, faces, landmarks(end))
        if i = =displ(k)
            options.start_points = landmarks
            subplot(2, 2, k)
            options.colorfx = 'equalize'
            plot_fast_marching_mesh(vertex, faces, min(Dland, [], 2) , [], options)
            k = k + 1


def exo2():
    """
    Create an interpolation scheme to interpolate the result of MDS
    dimensionality reduction with Stree minimization (SMACOF algorithm).
    """


