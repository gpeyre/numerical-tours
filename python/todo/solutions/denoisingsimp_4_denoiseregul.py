def exo1():
    """
    Compute the solution for several value of $\lambda$ and choose the
    optimal |lambda| and the corresponding optimal denoising |fSob0|.
    You can increase progressively lambda and reduce
    considerably the number of iterations.
    lot
    """
    lambda_list = linspace(1, 30, 60)
    err = []
    for i in 1: length(lambda_list):
        lambda = lambda_list(i)
        fSob = real(ifft2(yF ./ (1 + lambda*S)))
        err(i) = snr(f0, fSob)
    [tmp, i] = max(err)
    lambda = lambda_list(i)
    fSob0 = real(ifft2(yF ./ (1 + lambda*S)))
    plot(lambda_list, err); axis('tight')
    set_label('lambda', 'SNR')


def exo2():
    """
    Compute the total variation of |f0|.
    isplay
    """
    tv  = sum(sum(sqrt(sum(grad(f0).^2, 3))))
    imageplot(f0, strcat(['TV = ', num2str(tv, 4)]), 1, 2, 1)


def exo3():
    """
    Compute the gradient descent and monitor
    the minimized energy.
    """
    niter = 200
    fTV = y
    energy = []
    for i in 1: niter:
        Gr = grad(fTV)
        d = sqrt(sum3(Gr.^2, 3))
        deps = sqrt(epsilon^2 + d.^2)
        G0 = -div(Gr ./ repmat(deps, [1 1 2]))
        G = fTV-y + lambda*G0
        energy(i) = 1/ 2*norm(y-fTV, 'fro')^2 + lambda*sum(deps(: ))
        fTV = fTV - tau*G
    plot(1: niter, energy); axis('tight')
    set_label('iteration', 'Energy')


def exo4():
    """
    Compute the solution for several value of $\lambda$ and choose the
    optimal $\lambda$ and the corresponding optimal denoising |fSob0|.
    You can increase progressively $\lambda$ and reduce
    considerably the number of iterations.
    """
    niter = 800
    lambda_list = linspace(0, .25, niter)
    tau = 2 / (1 + max(lambda) * 8 / epsilon)
    fTV = y
    energy = []
    for i in 1: niter:
        lambda = lambda_list(i)
        Gr = grad(fTV)
        d = sqrt(sum3(Gr.^2, 3))
        G0 = -div(Gr ./ repmat(sqrt(epsilon^2 + d.^2) , [1 1 2]))
        G = fTV-y + lambda*G0
        deps = sqrt(epsilon^2 + d.^2)
        fTV = fTV - tau*G
        err(i) = snr(f0, fTV)
        if i >1
            if err(i) > max(err(1: i-1))
                fTV0 = fTV
    plot(lambda_list, err); axis('tight')
    set_label('\lambda', 'SNR')


def exo5():
    """
    Compare the TV denoising with a hard thresholding in a translation
    invariant tight frame of wavelets.
    """
    extend_stack_size(4)
    Jmin = 4
    options.ti = 1
    fW = perform_wavelet_transf(y, Jmin, + 1, options)
    fWT = perform_thresholding(fW, 2.9*sigma, 'hard')
    fWav = perform_wavelet_transf(fWT, Jmin, -1, options)
    ewav = snr(f0, fWav)
    imageplot(clamp(y), strcat(['Noisy ' num2str(enoisy, 3) 'dB']), 1, 2, 1)
    imageplot(clamp(fWav), strcat(['Wavelets TI ' num2str(ewav, 3) 'dB']), 1, 2, 2)


