def exo1():
    """
    Compute the ratio $M/N$ of non-zero coefficients.
    """
    M = sum(fWT(: )~ = 0)
    disp(['Ratio of non-zero coefficients: M/ N = ' num2str(M/ n0^2, 3) '.'])


def exo2():
    """
    Compute a threshold $T$ to keep only $M$ coefficients.
    """
    v = sort(abs(fW(: )))
    if v(1) <v(n0*n0)
        v = reverse(v)
    T = v(M + 1)


def exo3():
    """
    Compute an approximation with an decreasing number of coefficients.
    """
    M_list = round(n0^2 ./ [16 32 64 128])
    for i in 1: length(M_list):
        M = M_list(i)
        f1 = PsiS(Thresh(fW, v(M + 1)))
        imageplot(clamp(f1), strcat(['M/ N = 1/ ' num2str(n0^2/ M) ', SNR = ' num2str(snr(f, f1), 3) 'dB']), 2, 2, i)


def exo4():
    """
    Try to optimize the value of the threshold $T$ to get the best possible
    denoising result.
    """
    q = 30
    Tlist = linspace(2.5, 4, q)*sigma
    err = []
    for i in 1: q:
        fWT = Thresh(fW, Tlist(i))
        f1 = PsiS(fWT)
        err(i) = snr(f, f1)
    plot(Tlist/ sigma, err)
    axis('tight')
    set_label('T/ \sigma', 'SNR')


def exo5():
    """
    Compute the full denoising by cycling through the $i$ indices.
    """
    f1 = zeros(n0, n0)
    for i in 1: tau_max*tau_max:
        fTrans = circshift(y, [X(i) Y(i)])
        fTrans = PsiS(Thresh(Psi(fTrans) , T))
        fTrans = circshift(fTrans, -[X(i) Y(i)])
        f1 = (i-1)/ i*f1 + fTrans/ i
    imageplot(clamp(y), 'Noisy image', 1, 2, 1)
    imageplot(clamp(f1), strcat(['Denoising, SNR = ' num2str(snr(f, f1), 3) 'dB']), 1, 2, 2)


def exo6():
    """
    Determine the optimal threshold $T$ for this translation invariant
    denoising.
    """


def exo7():
    """
    Test on other images.
    """


