def exo1():
    """
    Implement the Dijkstra algorithm by iterating these step while the
    stack |I| is non empty.
    Display from time to time the front that propagates.
    """
    options.method = 'dijkstra'
    options.svg_rate = n*6
    [D, Dsvg, Ssvg] = perform_dijkstra_fm(W, x0, options)
    for i in 1: 4:
        subplot(2, 2, i)
        d = Dsvg(: , : , 2 + i); d(d = =Inf) = 0
        imageplot(d)
        colormap jet(256)


def exo2():
    """
    Implement the Fast Marching algorithm.
    Display from time to time the front that propagates.
    """
    options.method = 'fm'
    options.svg_rate = n*6
    [D, Dsvg, Ssvg] = perform_dijkstra_fm(W, x0, options)
    for i in 1: 4:
        subplot(2, 2, i)
        d = Dsvg(: , : , 2 + i); d(d = =Inf) = 0
        imageplot(d)
        colormap jet(256)


def exo3():
    """
    Compute the distance map to these starting point using the FM algorithm.
    """
    options.method = 'fm'
    D = perform_dijkstra_fm(W, x0, options)
    k = 8
    displ = lambda D: cos(2*pi*k*D/ max(D(: )))
    imageplot(displ(D))
    colormap jet(256)


def exo4():
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


