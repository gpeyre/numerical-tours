def exo1():
    """
    Display the surface with several light directions $d \in \RR^3$.
    """
    v = 2*rand(3, 4)-1; v(3, : ) = abs(v(3, : ))
    v = v ./ repmat(sqrt(sum(v.^2)), [3 1])
    for i in 1: 4:
        d = v(: , i)
        L = sum(N .* repmat(reshape(d, [1 1 3]), [n n 1]), 3)
        subplot(2, 2, i)
        imageplot(max(L, vmin))


def exo2():
    """
    Try to reconstruct the image starting from other base points.
    What do you observe ?
    """
    p = [[96; 72], [230; 175], [50; 125], [95; 180]]
    for i in 1: 1:
        [f1, S] = perform_fast_marching(1./ W, p(: , i))
        f1 = -f1*n
        subplot(1, 1, i)
        hold on
        surf(f1)
        h = plot3(p(2, i), p(1, i), f1(p(1, i), p(2, i)), 'r.')
        set(h, 'MarkerSize', 30)
        colormap(gray(256))
        shading interp
        axis('equal')
        view(110, 45)
        axis('off')
        camlight


def exo3():
    """
    Try to improve the quality of the reconstruction by selecting several
    points, and imposing their height.
    """


