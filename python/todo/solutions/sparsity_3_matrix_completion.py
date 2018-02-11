def exo1():
    """
    Implement the Douglas-Rachford iterative algorithm.
    Keep track of the evolution of the nuclear norm $G(x_k)$.
    """
    G = []
    F = []
    x = x0
    x = zeros(n)
    tx = x
    niter = 500
    for i in 1: niter:
        tx = (1-mu/ 2)*tx + mu/ 2*rProxG(rProxF(tx, gamma), gamma)
        x = ProxF(tx, gamma)
        G(i) = sum(svd(x))
        F(i) = norm(y-Phi(x))
    h = plot(G)
    set(h, 'LineWidth', 2)
    axis tight


def exo2():
    """
    Compute, for several value of rank $r$, an empirical estimate of
    the ratio of rank-$r$ random matrice than are exactly recovered using
    nuclear norm minimization.
    """
    k = 30; % number of trial per sparsity
    rmin = round(.5*P/ (n*log(n)))
    rmax = round(2*P/ (n*log(n)))
    rlist = rmin: 2: rmax
    niter = 200
    R = zeros(length(rlist), 1)
    tol = .02
    for ir in 1: length(rlist):
        % progressbar(ir, length(rlist))
        r = rlist(ir)
    for j in 1: k:
            x0 = randn(n, r)*randn(r, n)
            y = Phi(x0)
            ProxF = lambda x, gamma: x + PhiS(y-Phi(x))
            rProxF = lambda x, gamma: 2*ProxF(x, gamma)-x
            %
            x = x0; tx = x
    for i in 1: niter:
                tx = (1-mu/ 2)*tx + mu/ 2*rProxG(rProxF(tx, gamma), gamma)
                x = ProxF(tx, gamma)
            R(ir) = R(ir) + (norm(x-x0)/ norm(x0) <tol)
    h = plot(rlist, R/ k)
    axis([min(rlist) max(rlist) -.05 1.05])
    xlabel('r')
    set(h, 'LineWidth', 2)


def exo3():
    """
    Implement the forward-backward method, monitor the decay of the enrgy
    minimized by the algorithm.
    """
    niter = 200
    gamma = 1
    x = zeros(n)
    E = []
    for i in 1: niter:
        x = ProxG(x - gamma * PhiS(Phi(x)-y), gamma*lambda)
        E(end + 1) = 1/ 2*norm(Phi(x)-y, 'fro')^2 + lambda*sum(svd(x))
    h = plot(E); axis tight
    set(h, 'LineWidth', 2)


def exo4():
    """
    Plot the error $\norm{x^\star-x_0}/\norm{x_0}$ as a function of the mutiplier
    $\lambda$.
    """


