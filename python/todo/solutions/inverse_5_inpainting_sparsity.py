def exo1():
    """
    Perform the iterative soft thresholding.
    Monitor the decay of the energy $E$ you are minimizing.
    """
    fSpars = y
    energy = []
    niter = 1000
    for i in 1: niter:
        fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambda)
        % record the energy
        fW = PsiS(fSpars)
        energy(i) = 1/ 2 * norm(y-Phi(fSpars, Omega), 'fro')^2 + lambda * sum(abs(fW(: )))
    h = plot(energy)
    axis('tight')
    set_label('Iteration', 'E')
    if using_matlab()
        set(h, 'LineWidth', 2)


def exo2():
    """
    Since there is no noise, one should in theory takes $\lambda
    \rightarrow 0$.
    To do this, decay the value of $\lambda$ through the iterations.
    """
    niter = 1000
    lambda_list = linspace(.03, 0, niter)
    err = []
    for i in 1: niter:
        fSpars = SoftThreshPsi(ProjC(fSpars, Omega), lambda_list(i))
    imageplot(clamp(fSpars), ['Sparsity inpainting, SNR = ' num2str(snr(f0, fSpars), 3) 'dB'])


def exo3():
    """
    Perform the iterative soft thresholding. Monitor the decay of the
    energy $E$.
    """
    niter = 1000
    a = U.*PsiS(fSpars)
    E = []
    for i in 1: niter:
        fTI = Psi(a)
        d = y-Phi(fTI, Omega)
        E(i) = 1/ 2*norm(d , 'fro')^2 + lambda * sum(abs(a(: )))
        % step 
        a = SoftThresh(a + tau*PsiS(Phi(d, Omega)), lambda*tau)
    plot(E); axis('tight')


def exo4():
    """
    Perform the iteration with a decaying value of $\lambda$
    """
    niter = 3000
    lambda_list = linspace(.03, 0, niter)
    for i in 1: niter:
        fTI = Psi(a)
        d = y-Phi(fTI, Omega)
        % step 
        a = SoftThresh(a + tau * PsiS(Phi(d, Omega)) , lambda_list(i)*tau) 
    imageplot(clamp(fTI), ['Sparsity inpainting TI, SNR = ' num2str(snr(f0, fTI), 3) 'dB'])


def exo5():
    """
    Perform the iteration with a decaying value of $\lambda$
    """
    niter = 500
    lambda_list = linspace(1, 0, niter)
    fHard = y
    for i in 1: niter:
        fHard = Xi(HardThresh(PsiS(ProjC(fHard, Omega)), lambda_list(i)))
    imageplot(clamp(fHard), ['Inpainting hard thresh., SNR = ' num2str(snr(f0, fHard), 3) 'dB'])


