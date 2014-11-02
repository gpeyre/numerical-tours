def exo1():
    """
    Compute the denoising SNR for different values of |mu| and different
    value of |T|. Important: to get good results, you should not threshold
    the low frequency residual.
    list/sigma, mulist,
    """
    mulist = linspace(1, 12, 17)
    nmu = length(mulist)
    Tlist = linspace(.5, 3.5, 10)*sigma
    err = []
    for i in 1: length(mulist):
    for j in 1: length(Tlist):
            T = Tlist(j)
            MWT = perform_thresholding(MW, [T mulist(i)*T], 'semisoft')
            MWT(1: 2^Jmin, 1: 2^Jmin) = MW(1: 2^Jmin, 1: 2^Jmin)
            MT = perform_wavelet_transf(MWT, Jmin, -1)
            err(i, j) = snr(M0, MT)
    imageplot(err)
    axis('tight')
    set_label('T/ \sigma', '\mu')
    title('SNR of semi-soft thresholding')
    axis('on')


def exo2():
    """
    Compare the performance of Soft and Stein thresholders, by determining
    the best threshold value.
    """
    Tlist = linspace(.5, 4, 20)*sigma
    err_hard = []
    err_soft = []
    err_stein = []
    for j in 1: length(Tlist):
        T = Tlist(j)
        MWT = perform_thresholding(MW, T, 'hard')
        MWT(1: 2^Jmin, 1: 2^Jmin) = MW(1: 2^Jmin, 1: 2^Jmin)
        MT = perform_wavelet_transf(MWT, Jmin, -1)
        err_hard(j) = snr(M0, MT)
        MWT = perform_thresholding(MW, T, 'soft')
        MWT(1: 2^Jmin, 1: 2^Jmin) = MW(1: 2^Jmin, 1: 2^Jmin)
        MT = perform_wavelet_transf(MWT, Jmin, -1)
        err_soft(j) = snr(M0, MT)
        MWT = perform_thresholding(MW, T, 'stein')
        MWT(1: 2^Jmin, 1: 2^Jmin) = MW(1: 2^Jmin, 1: 2^Jmin)
        MT = perform_wavelet_transf(MWT, Jmin, -1)
        err_stein(j) = snr(M0, MT)
    hold('on')
    plot(Tlist/ sigma, err_hard, 'r')
    plot(Tlist/ sigma, err_soft, 'b')
    plot(Tlist/ sigma, err_stein, 'g')
    hold('off')
    axis('tight')
    set_label('T/ \sigma', 'SNR')
    legend('Hard', 'Soft', 'Stein')
    axis('on')


