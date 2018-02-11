def exo1():
    """
    Perform the iterative soft thresholding.
    Monitor the decay of the energy $E$ you are minimizing.
    """
    fSpars = y
    energy = []
    niter = 100
    for i in 1: niter:
        fSpars = fSpars + tau * Phi(y-Phi(fSpars))
        % thresholding
        fSpars = SoftThreshPsi(fSpars, lambda*tau)
        % record the energy
        fW = PsiS(fSpars)
        energy(end + 1) = 1/ 2 * norm(y-Phi(fSpars), 'fro')^2 + lambda * sum(abs(fW(: )))
    h = plot(energy); axis([1, niter, min(energy)*1.05 + max(energy)*(-.05) max(energy)])
    set_label('Iteration', 'E')
    set(h, 'LineWidth', 2)


def exo2():
    """
    Try to find the best threshold $\lambda$. To this end, perform a *lot*
    of iterations, and progressively decay the threshold $\lambda$ during the
    iterations. Record the best result in |fBestOrtho|.
    armup
    escent
    """
    switch setting
        case 1
            lambda_max = .01/ 4
        case 2
            lambda_max = .02
        otherwise
            error('Unknown setting #')
    niter = 1000
    fSpars = y
    lambda_list_Ortho = linspace(lambda_max, 0, niter)
    errOrtho = []
    for i in 1: niter:
        fSpars = fSpars + tau * Phi(y-Phi(fSpars))
        fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(1)*tau)
    for i in 1: niter:
        % descent
        fSpars = fSpars + tau * Phi(y-Phi(fSpars))
        % thresholding
        fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(i)*tau)
        % record the error
        errOrtho(i) = snr(f0, fSpars)
        if i >1 && errOrtho(i) >max(errOrtho(1: i-1))
            fBestOrtho = fSpars
    h = plot(lambda_list_Ortho, errOrtho)
    axis tight
    set_label('lambda', 'SNR')
    set(h, 'LineWidth', 2)


def exo3():
    """
    Perform the iterative soft thresholding. Monitor the decay of the
    energy.
    """
    niter = 500
    a = PsiS(y)
    E = []
    for i in 1: niter:
        fTI = Psi(a)
        d = y-Phi(fTI)
        % measure energy
        E(end + 1) = 1/ 2*norm(d, 'fro')^2 + lambda*sum(sum(sum(abs(a.*U))))
        % step
        a = SoftThresh(a + tau*PsiS(Phi(d)), lambda*tau)
    plot(E); axis('tight')


def exo4():
    """
    Compute the optimal value of $\lambda$, and record the optimal
    reconstruction |fBestTI|.
    armup
    escent
    """
    switch setting
        case 1
            lambda_max = .005
        case 2
            lambda_max = .01
        otherwise
            error('Unknown setting #')
    lambda_list_TI = linspace(lambda_max/ 4, lambda_max, niter)
    a = PsiS(y)*0
    errTI = []
    for i in 1: niter:
        fTI = Psi(a)
        d = y-Phi(fTI)
        a = SoftThresh(a + tau*PsiS(Phi(d)), lambda_list_TI(1)*tau)
    for i in 1: niter:
        fTI = Psi(a)
        d = y-Phi(fTI)
        % step 
        a = SoftThresh(a + tau*PsiS(Phi(d)), lambda_list_TI(i)*tau)
        errTI(end + 1) = snr(f0, fTI)
    	if i >1 && errTI(i) >max(errTI(1: i-1))
            fBestTI = fTI
    plot(lambda_list_TI, errTI); axis('tight')


def exo5():
    """
    Compare with the result of TV regularization, record the optimal
    TV result in |fBestTV|.
    armup stage
    escent
    lf; plot(log10(E(1:end/2)/E(end)-1));
    """
    switch setting
        case 1
            lambda_max = 0.002
        case 2
            lambda_max = .02
        otherwise
            error('Unknown setting #')
    epsilon = 1e-2
    niter = 1000
    lambda_list_TV = linspace(lambda_max/ 4, lambda_max, niter)
    tau = 1.9 / (1 + max(lambda_list_TV) * 8 / epsilon)
    fBestTV = y; fTV = y
    errTV = []; E = []
    for i in 1: 1000:
        lambda = lambda_list_TV(1)
        % Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV)
        d = sqrt(epsilon^2 + sum3(Gr.^2, 3))
        G = -div(Gr./ repmat(d, [1 1 2]))
        % step
        e = Phi(fTV)-y
        fTV = fTV - tau*(Phi(e) + lambda*G)
    for i in 1: niter:
        lambda = lambda_list_TV(i)
        % Compute the gradient of the smoothed TV functional.
        Gr = grad(fTV)
        d = sqrt(epsilon^2 + sum3(Gr.^2, 3))
        G = -div(Gr./ repmat(d, [1 1 2]))
        % step
        e = Phi(fTV)-y
        fTV = fTV - tau*(Phi(e) + lambda*G)
        % record error
        errTV(i) = snr(f0, fTV)
        E(i) = 1/ 2*norm(e, 'fro')^2 + lambda * sum(d(: ))
        if errTV(i) >snr(f0, fBestTV)
            fBestTV = fTV
    plot(lambda_list_TV, errTV)
    axis('tight')
    xlabel('\lambda'); ylabel('SNR')


