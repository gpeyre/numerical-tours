def exo1():
    """
    Compute the Lipschitz constant $L=\norm{\Ga_z}^2$. _Hint:_ you can use power iteration to estimate the largest eigenvalues
    of $\Ga_z^*\Ga_z$. We display here the evolution of the estimate of this eigenvalue with the number of power iterations.
    
    
    """
    A = lambda u: GammaS(Gamma(u))
    u = randn([N 2])
    u = u/ norm(u(: ))
    e = []
    for i in 1: 15:
        v = A(u)
        e(end + 1) = sum(u(: ).*v(: ))
        u = v/ norm(v(: ))
    L = e(end)
    clf; plot(e, 'LineWidth', 2); axis tight


def exo2():
    """
    Implement the forward-backward algorithm. Monitor the decay of the energy
    $ \log_{10}( E(u^{(\ell)})/ E^\star  -1   ) $
    where $E^\star$ is an estimate of the minimum energy (obtained by
    running the algorithm with a lot of iterations.
    """
    R = []
    u = zeros(N, 2)
    niter = 2000
    damp = 1.8
    tau = damp/ L
    for i in 1: niter:
        u = ProxJ(u - tau * gradF(u), lambda*tau)
        R(end + 1) = E(u)
    sel = 1: niter/ 4
    plot(sel, log(R(sel)/ min(R)-1), '-', 'LineWidth', 2); axis tight


def exo3():
    """
    Compute the full homotopy path $\la \mapsto (a_\la,b_\la) $
    where $ (a_\la,b_\la) $ is the solution of C-BP with regularization parameter
    $\la$. The display bellow shows in color the evolution of the correctly estimated spikes
    and in black wrong perturbating spikes.
    Test with different values of $\kappa$. What can you conclude about the performances of C-BP ?
    ist of regul param
    armup stage
    uter loop
    isplay
    """
    Lmax = max(PhiS(y))
    Llist = linspace(Lmax, 1e-6, 30)
    u = zeros(N, 2); % 
    niter = 3000
    lambda = Llist(1)
    for j in 1: niter:
    	u = ProxJ(u - tau * gradF(u), lambda*tau)
    niter = 500
    A = []; B = []
    for i in 1: length(Llist):
        % progressbar(i, length(Llist))
        lambda = Llist(i)
        % inner loop
    for j in 1: niter:
            u = ProxJ(u - tau * gradF(u), lambda*tau)
        a = u(: , 1);  b = u(: , 2)
        A(: , end + 1) = a; B(: , end + 1) = b
    Ic = setdiff(1: N, I)
    subplot(2, 1, 1); hold on
    plot(Llist, A(I, : )', 'LineWidth', lw)
    plot(Llist, A(Ic, : )', 'k-', 'LineWidth', lw); axis tight; box on
    title('a')
    subplot(2, 1, 2); hold on
    plot(Llist, B(I, : )', 'LineWidth', lw)
    plot(Llist, B(Ic, : )', 'k-', 'LineWidth', lw); axis tight; box on
    axis([0 max(Llist) -1 1])
    title('b')


