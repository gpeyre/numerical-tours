def exo1():
    """
    Compute and display the pseudo inverse reconstruction
    $ \Phi^+ y $. What do you observe ?
    """
    fL2 = PhiS(y)
    e = snr(f0, fL2)
    imageplot(clamp(fL2), ['SNR = ' num2str(e, 3) 'dB'])


def exo2():
    """
    Find the optimal solution |fSob| by testing several value of
    $\lambda$.
    """
    lambda_list = linspace(.01, 12, 40)
    err = []
    for i in 1: length(lambda_list):
        lambda = lambda_list(i)
        fSob = real(ifft2(y .* xi ./ (xi + lambda*S)))*n
        err(i) = snr(f0, fSob)
    h = plot(lambda_list, err); axis tight
    set_label('lambda', 'SNR')
    set(h, 'LineWidth', 2)
    [tmp, i] = max(err)
    lambda = lambda_list(i)
    fSob = real(ifft2(y .* xi ./ (xi + lambda*S)))*n


def exo3():
    """
    Perform the iterative soft thresholding.
    Monitor the decay of the energy $E$ you are minimizing.
    """
    fSpars = PhiS(y)
    energy = []
    niter = 100
    for i in 1: niter:
        fSpars = fSpars + tau * PhiS(y-Phi(fSpars))
        % thresholding
        fSpars = SoftThreshPsi(fSpars, lambda*tau)
        % record the energy
        fW = PsiS(fSpars)
        energy(end + 1) = 1/ 2 * norm(y-Phi(fSpars), 'fro')^2 + lambda * sum(abs(fW(: )))
    h = plot(energy); axis([1, niter, min(energy)*1.05 + max(energy)*(-.05) max(energy)])
    set_label('Iteration', 'E')
    set(h, 'LineWidth', 2)


def exo4():
    """
    Try to find the best threshold $\lambda$. To this end, perform a *lot*
    of iterations, and progressively decay the threshold $\lambda$ during the
    iterations. Record the best result in |fBestOrtho|.
    armup
    terations
    """
    lambda_max = .3
    niter = 200
    fSpars = PhiS(y)
    lambda_list_Ortho = linspace(lambda_max, 0, niter)
    errOrtho = []
    for i in 1: niter:
        fSpars = fSpars + tau * PhiS(y-Phi(fSpars))
        fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(1)*tau)
    for i in 1: niter:
        fSpars = fSpars + tau * PhiS(y-Phi(fSpars))
        fSpars = SoftThreshPsi(fSpars, lambda_list_Ortho(i)*tau)
        % record the error
        errOrtho(i) = snr(f0, fSpars)
        if i >1 && errOrtho(i) >max(errOrtho(1: i-1))
            fBestOrtho = fSpars
    h = plot(lambda_list_Ortho, errOrtho)
    axis tight
    set_label('lambda', 'SNR')
    set(h, 'LineWidth', 2)


def exo5():
    """
    Use the iterative thresholding but this time
    with the translation invariant wavelet transform.
    Find the best value of $\lambda$ and
    record the best result in |fBestTI|.
    armup
    terations
    """
    lambda_max = .3
    niter = 200
    fSpars = PhiS(y)
    lambda_list_TI = linspace(lambda_max, 0, niter)
    errTI = []
    for i in 1: niter:
        fSpars = fSpars + tau * PhiS(y-Phi(fSpars))
        fSpars = SoftThreshPsi(fSpars, lambda_list_TI(1)*tau)
    for i in 1: niter:
        fSpars = fSpars + tau * PhiS(y-Phi(fSpars))
        fSpars = SoftThreshPsi(fSpars, lambda_list_TI(i)*tau)
        % record the error
        errTI(i) = snr(f0, fSpars)
        if i >1 && errTI(i) >max(errTI(1: i-1))
            fBestTI = fSpars
    h = plot(lambda_list_TI, errTI)
    axis tight
    set_label('lambda', 'SNR')
    set(h, 'LineWidth', 2)


def exo6():
    """
    Compare Sobolev and Sparse reconstruction for MRI imaging.
    For a given number of Fourier sample, compare the quality of the
    reconstruction for different $\alpha$.
    """


