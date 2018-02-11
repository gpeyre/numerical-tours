def exo1():
    """
    Compute the full shape features |d(k,i)| for all points |i| and scales
    |k|.
    """
    [Y, X] = meshgrid(1: n, 1: n)
    d = []
    for k in 1: nrad:
        r = rlist(k)
        if 0
            Mh = perform_blurring(M1, r)
        else
            x = -ceil(r): ceil(r)
            [b, a] = meshgrid(x, x)
            h = double(a.^2 + b.^2 <= r^2)
            h = h/ sum(h(: ))
            Mh = perform_convolution(M, h)
        [Y, X] = meshgrid(1: n, 1: n)
    for i in 1: 2:
            D{i}(k, : ) = interp2(Y, X, Mh{i}, bound{i}(2, : ), bound{i}(1, : ))
            % I = round(bound{i}(1, : )) + round(bound{i}(2, : )-1)*n
            % D{i}(k, : ) = Mh{i}(I)
    sel = 1: round(nbound/ 8)
    plot(sel, D{1}([1 nrad/ 2 nrad], sel)', '-')
    title('Some features (zoom)')
    axis tight


def exo2():
    """
    Compute the descriptor for all the values of |r| in rlist.
    """
    C = zeros(nbound)
    for i in 1: nbound:
    for j in 1: nbound:
            C(i, j) = norm(D{1}(: , i)-D{2}(: , j))
    imageplot(C)


def exo3():
    """
    Compute a metric associated to |C|, by rescaling.
    Using the Fast Marching, compute the shortest paths |gpath| from point |(1,1)|
    to point |(nbound,nbound)|. Record the length of this path, which is
    the value of the geodesic distance.
    etric
    istance
    ath
    """
    W = rescale(C, .001, 1)
    [D, S] = perform_fast_marching(1./ W, [1; 1])
    gpath = compute_geodesic(D, size(D))


