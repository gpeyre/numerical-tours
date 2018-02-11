def exo1():
    """
    Iterate many time the randomized assignement until convergence of $\tilde f$.
    The random projector $\Theta$ should be re-computed at each iteration.
    """
    niter = 1000
    f1 = f
    disp_list = [1 3 10 100]; q = 1
    for i in 1: niter:
        [Theta, ~] = qr(randn(d))
        f1 = (1-tau)*f1 + tau * Theta * P(Theta'*f1, Theta'*g)
        if q <= 4 && i = =disp_list(q)
            t = (q-1)/ 3
            subplot(2, 2, q)
            hold on
            plotp(f1, [t 0 1-t])
            axis('off'); axis('equal')
            q = q + 1


def exo2():
    """
    Show the progressive interpolation for varying $t \in [0,1]$.
    """
    tlist = linspace(0, 1, 6)
    for i in 1: length(tlist):
        t = tlist(i)
        ft = (1-t)*f + t*f1
        subplot(2, length(tlist)/ 2, i)
        plotp(ft, [t 0 1-t])
        axis('off'); axis('equal')
        title(['t = ' num2str(t, 2)])


def exo3():
    """
    Perform the equalization of each of the coordinate independantly
    of $f$ with $g$. Display the resulting image.
    """
    f1 = P(f, g)
    F1 = reshape(f1', [n n 3])
    imageplot(F1)


def exo4():
    """
    To obtain an exact matching, one can use the stochastic gradient descent
    algorithm to minimize $SW$.
    Display the resulting image at several stages of the optimization
    process.
    
    """
    niter = 1000
    f1 = f
    q = 1; disp_list = [3 10 100 niter]
    for i in 1: niter:
        [Theta, ~] = qr(randn(d))
        f1 = (1-tau)*f1 + tau * Theta * P(Theta'*f1, Theta'*g)
        if q <= 4 && i = =disp_list(q)
            subplot(2, 2, q)
            F1 = reshape(f1', [n n 3])
            imageplot(F1)
            q = q + 1
    F1 = reshape(f1', [n n 3])


def exo5():
    """
    Display the geodesic interpolation between the two histograms $\mu_f$ and $\mu_g$.
    """
    tlist = linspace(0, 1, 6)
    for i in 1: length(tlist):
        t = tlist(i)
        ft = (1-t)*f + t*f1
        subplot(2, length(tlist)/ 2, i)
        imageplot(func(hist2d(ft(1: 2, : ), Q)))
        title(['t = ' num2str(t, 2)])


