def exo1():
    """
    Compute the gradient descent and monitor
    the minimized energy.
    """
    niter = 200
    u = f
    energy = []
    for i in 1: niter:
        energy(i) = 1/ 2*norm(f-u, 'fro')^2 + lambda*J(u)
        u = u - tau*(u - f + lambda*GradJ(u))
    plot(1: niter, energy); axis('tight')
    set_label('iteration', 'Energy')


def exo2():
    """
    Perform the gradient descent, monitor the decay of the energy.
    """
    niter = 500
    u = f
    E = []
    for i in 1: niter:
        u = u - tau * (GradT(u-f) + lambda*GradJ(u))
        E(end + 1) = T(f-u) + lambda*J(u)
    plot(E)
    axis('tight')


def exo3():
    """
    Define the operators $\text{Grad} T$
    and apply it to an images.
    """
    GradT = lambda f: real(ifft2(W.^2.*fft2(f)))
    T = lambda v: 1/ 2*norm(W.*fft2(v)/ n, 'fro').^2
    imageplot(GradT(f), 'Grad(T)', 1, 2, 1)
    imageplot(f-GradT(f), 'f-Grad(T)', 1, 2, 2)


def exo4():
    """
    For a well chosen value of $\lambda$, perform the TV-Hilbert
    decomposition with this texture kernel.
    
    """
    lambda = .05
    tau = 1.9 / (max(W(: ))^2 + 8*lambda/ epsilon)
    niter = 500
    E = []
    u = f
    for i in 1: niter:
        u = u - tau * (GradT(u-f) + lambda*GradJ(u))
        E(end + 1) = T(f-u) + lambda*J(u)
    plot(E)
    axis('tight')


