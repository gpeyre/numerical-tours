def exo1():
    """
    Compute the solution to the heat equation.
    """
    f = f0
    clf; k = 0
    for i in 1: niter:
        f = f + tau * delta(f)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            imageplot(clamp(f), strcat(['T = ' num2str(T*k/ 4, 3)]), 2, 2, k)


def exo2():
    """
    Display the heat convolution for increasing values of $t$.
    """
    tlist = linspace(0.5, 10, 4)
    for i in 1: length(tlist):
        t = tlist(i)
        imageplot(heat(f0, t), ['t = ' num2str(t, 2)], 2, 2, i)


