def exo1():
    """
    Determine the best threshold $T$ for both hard and soft thresholding.
    Test several $T$ values in $[.8*\sigma, 4.5\sigma$,
    and display the empirical SNR $-10\log_{10}(\norm{f_0-\tilde f}/\norm{f_0})$
    What can you conclude from these results ?
    Test with another image.
    """
    Tlist = linspace(.8, 4.5, 25)*sigma
    err_soft = []; err_hard = []
    for i in 1: length(Tlist):
        aT = perform_thresholding(a, Tlist(i), 'hard')
    	fWav = perform_wavelet_transf(aT, Jmin, -1, options)
        err_hard(i) = snr(f0, fWav)
        aT = perform_thresholding(a, Tlist(i), 'soft')
        aT(1: 2^Jmin, 1: 2^Jmin) = a(1: 2^Jmin, 1: 2^Jmin)
    	fWav = perform_wavelet_transf(aT, Jmin, -1, options)
        err_soft(i) = snr(f0, fWav)
    h = plot(Tlist/ sigma, [err_hard(: ) err_soft(: )]); axis('tight')
    if using_matlab()
        set(h, 'LineWidth', 2)
    set_label('T/ \sigma', 'SNR')
    legend('Hard', 'Soft')


def exo2():
    """
    Perform the cycle spinning denoising by iterating on $i$.
    """
    fTI = zeros(n, n)
    T = 3*sigma
    for i in 1: m^2:
        fS = circshift(f, [dX(i) dY(i)])
        a = perform_wavelet_transf(fS, Jmin, 1, options)
        aT = perform_thresholding(a, T, 'hard')
        fS = perform_wavelet_transf(aT, Jmin, -1, options)
        fS = circshift(fS, -[dX(i) dY(i)])
        fTI = (i-1)/ i*fTI + 1/ i*fS
    imageplot(clamp(fHard), strcat(['Hard denoising, SNR = ' num2str(snr(f0, fHard), 3)]), 1, 2, 1)
    imageplot(clamp(fTI), strcat(['Cycle spinning denoising, SNR = ' num2str(snr(f0, fTI), 3)]), 1, 2, 2)


def exo3():
    """
    Study the influence of the number $m$ of shift on the denoising
    quality.
    """
    err = []
    shift_list = 1: 7
    T = 3*sigma
    for m in shift_list:
        [dY, dX] = meshgrid(0: m-1, 0: m-1)
        fTI = zeros(n, n)
    for i in 1: m^2:
            fS = circshift(f, [dX(i) dY(i)])
            a = perform_wavelet_transf(fS, Jmin, 1, options)
            aT = perform_thresholding(a, T, 'hard')
            fS = perform_wavelet_transf(aT, Jmin, -1, options)
            fS = circshift(fS, -[dX(i) dY(i)])
            fTI = (i-1)/ i*fTI + 1/ i*fS
        err(m) = snr(f0, fTI)
    h = plot(shift_list, err, '.-')
    if using_matlab()
        set(h, 'LineWidth', 2)
    axis('tight')
    set_label('m', 'SNR')


def exo4():
    """
    Determine the best threshold $T$ for both hard and soft thresholding,
    but now in the translation invariant case. What can you conclude ?
    """
    Tlist = linspace(.8, 4.5, 15)*sigma
    err_soft = []; err_hard = []
    for i in 1: length(Tlist):
        aT = perform_thresholding(a, Tlist(i), 'hard')
    	fWav = perform_wavelet_transf(aT, Jmin, -1, options)
        err_hard(i) = snr(f0, fWav)
        aT = perform_thresholding(a, Tlist(i), 'soft')
        aT(1: 2^Jmin, 1: 2^Jmin) = a(1: 2^Jmin, 1: 2^Jmin)
    	fWav = perform_wavelet_transf(aT, Jmin, -1, options)
        err_soft(i) = snr(f0, fWav)
    plot(Tlist(: )/ sigma, [err_hard(: ) err_soft(: )]); axis('tight')
    set_label('T/ \sigma', 'SNR')
    legend('Hard', 'Soft')


