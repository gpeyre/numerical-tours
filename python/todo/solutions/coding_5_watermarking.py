def exo1():
    """
    Compute |rho| so that |PSNR(y,x0,1)=snr_embedding|.
    """
    rho = sqrt(P)/ norm(w.*abs(x0)) * 10^(-psnr_embedding/ 20)
    disp(['rho = ' num2str(rho, 3) '.'])


def exo2():
    """
    According to you, for which PSNR the watermark becomes unoticeable?
    """


def exo3():
    """
    Using a Monte Carlo simulation (generation of the order of $10^3$
    watermarks, display the histogram of the repartition of $C(x_0,w)$.
    Compute the variance $\sigma_0^2$ of this distribution.
    """
    c = []
    q = 100
    niter = 60
    for i in 1: niter:
        c = [c C(repmat(x0, [1 q]), randn(P, q))]
    sigma0 = std(c)
    t = linspace(-4*sigma0, 4*sigma0, 31)
    h = hist(c, t)
    bar(t, h); axis('tight')


def exo4():
    """
    Compare, for various values of $T$ the estimation obtained by
    the Gaussian approximation with the true value obtained with the
    incomplete beta function.
    """


def exo5():
    """
    Compute, by Monte Carlo sampling (i.e. draw at random many $w$)
    the distribution of $C(A(x),w)$ for $x = x_0 + \rho \abs{x_0} w$. Store the different realization of
    $C(A(x),w)$ in a vector |c|.
    _Note:_ the value of $\rho$ should
    be recomputed for each $w$.
    """
    niter = 4000
    c = []
    for i in 1: niter:
        w = randn(P, 1)
        rho = sqrt(P)/ norm(w.*abs(x0)) * 10^(-psnr_embedding/ 20)
        x = x0 + rho*abs(x0).*w
        c(end + 1) = C(A(x), w)
    hist(c, 30)


def exo6():
    """
    Compute, for a varying value of $ p_{\text{FA}} $, the corresponding
    value of $ p_{\text{TP}} $. Display the resulting curve (ROC curve).
    This computation should be performed experimentally
    using e.g. 1000 random sampling.
    """
    pfa = 10.^(linspace(-10, -2))
    ptp = []
    for i in 1: length(pfa):
        T = sqrt(2) * sigma0 * erfinv(1-2*pfa(i))
        ptp(i) = sum(c >T)/ length(c)
    h = plot(log10(pfa), ptp)
    set(h, 'LineWidth', 2)
    axis('tight')
    xlabel('log_{10}(p_{FA})')
    ylabel('p_{TP}')


def exo7():
    """
    Try different attack strengths, by changing the value of $\tau$.
    For a $p_{\text{FA}}=10^{-6}$, determine the value of $\tau$
    for witch $p_{\text{TP}}$ drops bellow $0.2$.
    """


def exo8():
    """
    Try different attacks, for instance on the image itself (blurring,
    denoising, etc.).
    """


