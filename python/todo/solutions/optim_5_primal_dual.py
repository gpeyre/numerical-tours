def exo0():
    """
    Implement the primal-dual algorithm.
    Monitor the evolution of the TV energy $F(K(f_k))$
    during the iterations.
    Note that one always has $ f_k \in H $ so that the iterates
    satisfies the constraints.
    """
    niter = 200
    E = []; C = []
    for i in 1: niter:
        % update
        fold = f
        g = ProxFS(g + sigma*K(f1), sigma)
        f = ProxG(f-tau*KS(g), tau)
        f1 = f + theta * (f-fold)
        % monitor the decay of the energy
        E(i) = F(K(f))
        C(i) = snr(f0, f)
    h = plot(E)
    set(h, 'LineWidth', 2)
    axis('tight')


def exo1():
    """
    Use the primal dual scheme to perform regularization in the presence of
    noise
    $$ \umin{\norm{y-\Phi(f)} \leq \epsilon} \norm{\nabla f}_1. $$
    """


def exo2():
    """
    Display the evolution of the inpainting process.
    """
    y = Phi(f0)
    ProxG = lambda f, tau: f + Phi(y - Phi(f))
    niter = 600
    ndisp = round(linspace(1, niter, 5)); ndisp(1) = []
    E = []
    f = y
    g = K(y)*0
    f1 = f
    q = 1
    for i in 1: niter:
        % update
        fold = f
        g = ProxFS(g + sigma*K(f1), sigma)
        f = ProxG(f-tau*KS(g), tau)
        f1 = f + theta * (f-fold)
        % monitor the decay of the energy
        E(i) = F(K(f))
        if i = =ndisp(q)
            subplot(2, 2, q)
            imageplot(f)
            q = q + 1


