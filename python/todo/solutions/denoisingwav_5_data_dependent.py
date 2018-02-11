def exo1():
    """
    Display noisy image contaminated by Poisson noise of varying range.
    """
    lmin = 1
    lmax = [5 10 50 100]
    for i in 1: length(lmax):
        f1 = poissrnd(floor(rescale(f0u, lmin, lmax(i))))
        imageplot(f1, strcat(['\lambda_{max} = ' num2str(lmax(i))]), 2, 2, i)


def exo2():
    """
    Perform translation invariant wavelet hard thresholding directly on the
    Poisson noisy image $f$. Check for an optimal threshold that maximize
    the SNR. Record the optimal result in |fPoisson|.
    isplay
    est result
    """
    ntest = 30
    sigma = std(f(: )-f0(: ))
    Tlist = linspace(2, 3.5, ntest)*sigma
    fW = perform_wavelet_transf(f, Jmin, + 1, options)
    err = []
    for i in 1: ntest:
        fWT = perform_thresholding(fW, Tlist(i), 'hard')
        f1 = perform_wavelet_transf(fWT, Jmin, -1, options)
        err(i) = snr(f0, f1)
    h = plot(Tlist/ sigma, err)
    axis('tight')
    if using_matlab()
        set(h, 'LineWidth', 2)
    set_label('T/ sigma', 'SNR')
    [tmp, i] = max(err)
    fWT = perform_thresholding(fW, Tlist(i), 'hard')
    fPoisson = perform_wavelet_transf(fWT, Jmin, -1, options)


def exo3():
    """
    Display the estimated variance of a Poisson distribution for various
    $\lambda$ (e.g. 10000 realization for each $\lambda$) and display
    the variance of a stabilized distribution (here the green
    curve corresponds to 'Freeman & Tukey' and the blue curve to 'Anscombe'.
    lot
    bplot(2,1,1);
    = plot(v); axis('tight');
    using_matlab()
     set(h, 'LineWidth', 2);
    d
    tle('Poisson variance');
    t_label('\lambda', 'Variance');
    bplot(2,1,2);
    """
    lmax = 10
    ntest = 100000
    [V, U] = meshgrid(1: lmax, ones(ntest, 1))
    W = poissrnd(V)
    Mstab1 = 2*sqrt(W + 3/ 8)
    Mstab2 = sqrt(W + 1) + sqrt(W)
    vstab = []
    vstab(1, : ) = std(Mstab1, 1).^2
    vstab(2, : ) = std(Mstab2, 1).^2
    v = std(W, 1).^2
    h = plot(vstab'); axis('tight')
    if using_matlab()
        set(h, 'LineWidth', 2)
        legend('Anscombe', 'Freeman & Tukey')
    set_label('\lambda', 'Variance')


def exo4():
    """
    Perform translation invariance wavelet hard thresholding on the
    variance stabilized image. Use for instance the Anscombe VST.
    Check for an optimal threshold that maximize
    the SNR. Record the optimal result in |fVST|.
    isplay
    est result
    """
    ntest = 30
    sigma = 1
    Tlist = linspace(2, 3.5, ntest)*sigma
    fW = perform_wavelet_transf(2*sqrt(f + 3/ 8) , Jmin, + 1, options)
    err = []
    for i in 1: ntest:
        fWT = perform_thresholding(fW, Tlist(i), 'hard')
        f1 = perform_wavelet_transf(fWT, Jmin, -1, options)
        % undo VST
        f1 = (f1/ 2).^2 - 3/ 8
        % record error
        err(i) = snr(f0, f1)
    h = plot(Tlist/ sigma, err)
    axis('tight')
    if using_matlab()
        set(h, 'LineWidth', 2)
    set_label('T/ sigma', 'SNR')
    [tmp, i] = max(err)
    fWT = perform_thresholding(fW, Tlist(i), 'hard')
    f1 = perform_wavelet_transf(fWT, Jmin, -1, options)
    fVST = (f1/ 2).^2 - 3/ 8


def exo5():
    """
    Perform VST denoising using the Freeman VST.
    """


def exo6():
    """
    Generate several noisy images for several noise levels.
    """
    slist = [.1 .2 .3 .6]
    for i in 1: length(slist):
        Wu = gamrnd(1/ slist(i)^2, slist(i)^2, n, n)
        imageplot(f0.*Wu, strcat(['\sigma = ' num2str(slist(i))]), 2, 2, i)


def exo7():
    """
    Perform translation invariance wavelet hard thresholding directly on the
    noisy image $f=f_0 W$. Check for an optimal threshold that maximize
    the SNR. Record the optimal result in |fMult|.
    isplay
    est result
    """
    ntest = 30
    sigma = std(f(: )-f0(: ))
    Tlist = linspace(2.8, 4.5, ntest)*sigma
    fW = perform_wavelet_transf(f, Jmin, + 1, options)
    err = []
    for i in 1: ntest:
        fWT = perform_thresholding(fW, Tlist(i), 'hard')
        f1 = perform_wavelet_transf(fWT, Jmin, -1, options)
        err(i) = snr(f0, f1)
    h = plot(Tlist/ sigma, err)
    axis('tight')
    if using_matlab()
        set(h, 'LineWidth', 2)
    set_label('T/ sigma', 'SNR')
    [tmp, i] = max(err)
    fWT = perform_thresholding(fW, Tlist(i), 'hard')
    fMult = perform_wavelet_transf(fWT, Jmin, -1, options)


def exo8():
    """
    Perform translation invariance wavelet hard thresholding on the
    variance stabilized image using the log.
    Check for an optimal threshold that maximize
    the SNR. Record the optimal result in |fVST|.
    isplay
    est result
    """
    ntest = 30
    sigma = std(log(W(: )))
    Tlist = linspace(2, 3.5, ntest)*sigma
    fW = perform_wavelet_transf(log(f)-a, Jmin, + 1, options)
    err = []
    for i in 1: ntest:
        fWT = perform_thresholding(fW, Tlist(i), 'hard')
        f1 = perform_wavelet_transf(fWT, Jmin, -1, options)
        % undo VST
        f1 = exp(f1)
        % record error
        err(i) = snr(f0, f1)
    h = plot(Tlist/ sigma, err)
    axis('tight')
    if using_matlab()
        set(h, 'LineWidth', 2)
    set_label('T/ sigma', 'SNR')
    [tmp, i] = max(err)
    fWT = perform_thresholding(fW, Tlist(i), 'hard')
    f1 = perform_wavelet_transf(fWT, Jmin, -1, options)
    fVST = exp(f1)


