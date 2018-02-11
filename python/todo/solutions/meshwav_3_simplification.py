def exo1():
    """
    Perform iterative collapse to reach |p = round(2*n/3)| vertices.
    isplay
    """
    p = round(n/ 3)
    faces1 = faces
    vertex1 = vertex
    for i in 1: n-p:
        edges = compute_edges(faces1)
        nedges = size(edges, 2)
        k = floor(rand*(nedges-1)) + 1
        e = edges(: , k)
        vertex1(: , e(1)) = mean(vertex1(: , e), 2)
        vertex1(: , e(2)) = Inf
        faces1(faces1 = =e(2)) = e(1)
        a = sum(diff(sort(faces1)) = =0)
        faces1(: , a >0) = []
    plot_mesh(vertex1, faces1, options)
    shading faceted


def exo2():
    """
    As a post processing, find a way to remove from |faces1| and |vertex1| the unecessary
    information (remove vertex and faces that are not used).
    """


def exo3():
    """
    Perform iterative collapse to reach |p = round(2*n/3)| vertices.
    Use an ordering of the edge by their length.
    isplay
    """
    p = round(n/ 3)
    faces1 = faces
    vertex1 = vertex
    for i in 1: n-p:
        edges = compute_edges(faces1)
        D = vertex(: , edges(1, : )) - vertex(: , edges(2, : ))
        D = sum(D.^2, 1)
        [tmp, k] = min(D)
        e = edges(: , k)
        vertex1(: , e(1)) = mean(vertex1(: , e), 2)
        vertex1(: , e(2)) = Inf
        faces1(faces1 = =e(2)) = e(1)
        a = sum(diff(sort(faces1)) = =0)
        faces1(: , a >0) = []
    plot_mesh(vertex1, faces1, options)
    shading faceted


def exo4():
    """
    Try to use other criteria.
    o correction for this exercise.
    """


