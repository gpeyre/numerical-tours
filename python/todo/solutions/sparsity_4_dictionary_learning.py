def exo1():
    """
    Perform the iterative hard thresholding,
    and display the decay of the energy $J(x_j) = \norm{y_j-D x_j}^2$ for several $j$.
    _Remark:_ note that the iteration can be performed in parallel on all
    $x_j$.
    """
    niter = 100
    gamma = 1.6/ norm(D)^2
    E = []
    X = zeros(p, m)
    for i in 1: niter:
        R = D*X-Y
        E(end + 1, : ) = sum(R.^2)
        X = ProjX(X - gamma * D'*R, k)
    sel = 1: 5
    plot(log10(E(1: end/ 2, sel) - repmat(min(E(: , sel), [], 1), [niter/ 2 1])))
    axis tight
    title('log_{10}(J(x_j) - J(x_j^*))')


def exo2():
    """
    Perform this gradient descent, and monitor the decay of the energy.
    """
    niter_dico = 100
    E = []
    tau = 1/ norm(X*X')
    for i in 1: niter_dico:
        R = D*X - Y
        E(end + 1) = sum(R(: ).^2)
        D = ProjC(D + tau * (Y-D*X)*X')
    plot(log10(E(1: end/ 2)-min(E)))
    axis tight


def exo3():
    """
    Perform the dictionary learning by iterating between sparse coding and
    dictionary update.
    """
    niter_learning = 10
    niter_dico = 50
    niter_coef = 100
    E0 = []
    X = zeros(p, m)
    D = D0
    for iter  in  1: niter_learning:
        % --- coefficient update ----
        E = []
        gamma = 1.6/ norm(D)^2
    for i in 1: niter:
            R = D*X - Y
            E(end + 1, : ) = sum(R.^2)
            X = ProjX(X - gamma * D'*R, k)
        E0(end + 1) = norm(Y-D*X, 'fro')^2
        % --- dictionary update ----
        E = []
        tau = 1/ norm(X*X')
    for i in 1: niter_dico:
            R = D*X - Y
            E(end + 1) = sum(R(: ).^2)
            D = ProjC(D - tau * (D*X - Y)*X')
        E0(end + 1) = norm(Y-D*X, 'fro')^2
    clf; hold on
    plot(1: 2*niter_learning, E0)
    plot(1: 2: 2*niter_learning, E0(1: 2: 2*niter_learning), '*')
    plot(2: 2: 2*niter_learning, E0(2: 2: 2*niter_learning), 'o')
    axis tight
    legend('|Y-DX|^2', 'After coefficient update', 'After dictionary update')


