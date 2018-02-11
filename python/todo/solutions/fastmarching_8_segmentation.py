def exo1():
    """
    Display geodesic balls {x \ M(x)<T} for various T.
    """
    v = sort(D(: ))
    Tlist = v(round([.05 .1 .15 .25]*n^2))
    for i in 1: 4:
        T = Tlist(i)
        A = repmat(M, [1 1 3])
        I = find(D <T)
        A(I) = 1
        A([I + n^2; I + 2*n^2]) = 0
        subplot(2, 2, i)
        hold on
        imageplot(A)
        [c, h] = contour(D <T, [.5 .5], 'b')
        set(h, 'LineWidth', 2.5)
        axis('ij')


def exo2():
    """
    Display the level sets.
    """
    col = {'r' 'g' 'b' 'c' 'y' 'k' 'r: '}
    hold on
    imageplot(M)
    Vlist = unique(Q(: ))'
    for i in Vlist:
        U = zeros(n); U(Q = =i) = 1
        [c, h] = contour(U, [.5 .5], col{i})
        set(h, 'LineWidth', 4)
    axis ij
    hold off


