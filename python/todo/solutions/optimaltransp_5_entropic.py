def exo1():
    """
    Perform the iterations, and display the decay of the errors
    $$ \norm{\pi_\ell \ones - p}
     \qandq
     \norm{\pi_\ell^* \ones - q} $$
    in log scale.
    isplay error decay.
    """
    u = ones(N, 1)
    niter = 500
    gamma = .001
    pi = exp(-C/ gamma)
    E1 = []; E2 = []
    for i in 1: niter:
        pi = ProjC1(pi, p)
        E2(end + 1) = norm(pi'*u-q)/ norm(q)
        pi = ProjC2(pi, q)
        E1(end + 1) = norm(pi*u-p)/ norm(p)
    plot(log10([E1; E2]'))
    axis tight


def exo2():
    """
    Display the transport map for several values of $\gamma$.
    """
    glist = [.1 .01 .001 .0001]
    niter = 500
    for ig in 1: length(glist):
        gamma = glist(ig)
        pi = exp(-C/ gamma)
    for i in 1: niter:
            pi = ProjC2(ProjC1(pi, p), q)
        imageplot(normalizeMax(pi), ['\gamma = ' num2str(gamma)], 2, 2, ig)


