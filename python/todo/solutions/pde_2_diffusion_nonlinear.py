def exo1():
    """
    Implement the Perona-Malick diffusion flow
    for $\la = 10^{-2}$.
    """
    lambda = 1e-2
    T = .5/ lambda
    tau = .2
    niter = ceil(T/ tau)
    f = f0
    clf; k = 0
    for i in 1: niter:
        f = f + tau * div(g(A(f), lambda) .* grad(f))
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            imageplot(clamp(f), strcat(['T = ' num2str(T*k/ 4, 3)]), 2, 2, k)


def exo2():
    """
    Implement the Perona-Malick diffusion flow
    for $\la = 10^{-3}$.
    """
    lambda = 1e-3
    T = .5/ lambda
    tau = .5
    niter = ceil(T/ tau)
    f = f0
    clf; k = 0
    for i in 1: niter:
        f = f + tau * div(g(A(f), lambda) .* grad(f))
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            imageplot(clamp(f), strcat(['T = ' num2str(T*k/ 4, 3)]), 2, 2, k)


def exo3():
    """
    Implement the mean curvature flow.
    """
    T = 100
    tau = .2
    niter = ceil(T/ tau)
    f = f0
    clf; k = 0
    for i in 1: niter:
        f = f + tau * amplitude(grad(f)) .* curv(f)
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            imageplot(clamp(f), strcat(['T = ' num2str(T*k/ 4, 3)]), 2, 2, k)


def exo4():
    """
    Implement the affine-invariant curvature flow.
    """
    E = lambda s: sign(s).*abs(s).^(1/ 3)
    T = 100
    tau = .2
    niter = ceil(T/ tau)
    f = f0
    clf; k = 0
    for i in 1: niter:
        f = f + tau * amplitude(grad(f)) .* E(curv(f))
        if mod(i, floor(niter/ 4)) = =0
            k = k + 1
            imageplot(clamp(f), strcat(['T = ' num2str(T*k/ 4, 3)]), 2, 2, k)


