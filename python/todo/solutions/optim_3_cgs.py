def exo1():
    """
    Implement the conjugate gradient method, and monitor the decay of the
    energy $\norm{r_k}=\norm{Ax_k-b}$.
    """
    niter = 300
    x = zeros(n, 1)
    r = b-A*x
    p = r
    E = []
    for i in 1: niter:
        r1 = r
        E(i) = norm(r(: ))
        Ap = A*p
        alpha = dotp(r, r) / dotp(p, Ap)
        x = x + alpha*p
        r = r-alpha*Ap
        beta = dotp(r, r)/ dotp(r1, r1)
        p = r + beta*p
    plot(log10(E))
    title('log_{10}|r_k|')
    axis tight


def exo2():
    """
    Implement the conjugate gradient method, and monitor the decay of the
    energy $F(x_k)$.
    """
    niter = 200
    r = b
    p = r
    x = zeros(n)
    E = []; F = []
    for i in 1: niter:
        r1 = r
        E(i) = norm(r(: ))
        Ap = A(p)
        alpha = dotp(r, r) / dotp(p, Ap)
        x = x + alpha*p
        r = r-alpha*Ap
        beta = dotp(r, r)/ dotp(r1, r1)
        p = r + beta*p
    plot(log10(E))
    title('log_{10}|r_k|')
    axis tight


def exo3():
    """
    Implement the conjugate gradient method, and monitor the decay of the
    energy $F(x_k)) = \norm{\nabla x_k}$
    and the constraint $C(x_k) = \norm{y-\Phi x_k}^2$.
    _Important:_ be carefull at the initialization of the method.
    """
    niter = 600
    z = randn(n, n, 2)
    r = b-A(z)
    p = r
    C = []; F = []
    for i in 1: niter:
        %
        x = z(: , : , 1)
        C(i) = norm(y-Phi(x))
        F(i) = norm(delta(x))
        %
        r1 = r
        Ap = A(p)
        alpha = dotp(r, r) / dotp(p, Ap)
        z = z + alpha*p
        r = r-alpha*Ap
        beta = dotp(r, r)/ dotp(r1, r1)
        p = r + beta*p
    subplot(2, 1, 1)
    plot(log10(C))
    title('log_{10}|C(x_k)|')
    axis tight
    subplot(2, 1, 2)
    plot(log(abs(F(1: end/ 2)-F(end))))
    title('log_{10}|F(x_k) - F(x^{*})|')
    axis tight


