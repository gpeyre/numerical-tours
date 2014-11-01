    """
    Compute the distance map to these starting point using the FM algorithm.
    """
    options.method = 'fm'
    D = perform_dijkstra_fm(W, x0, options)
    k = 8
    displ = lambda D: cos(2*pi*k*D/ max(D.flatten()))
    imageplot(displ(D))
    colormap(jet(256))