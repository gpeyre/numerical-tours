def exo1():
    """
    Display the evolution of the denoising SNR when $T$ varies.
    Store in |fThresh| the optimal denoising result.
    etrieve best
    """
    Tlist = linspace(.5, 2.5, 30)*sigma
    snr_stein = arrayfun(lambda T: snr(f0, ThreshWav(f, T)), Tlist)
    clf; plot(Tlist/ sigma, snr_stein(: ), 'LineWidth', 2)
    axis tight; set_label('T/ \sigma', 'SNR')
    [~, i] = max(snr_stein); T = Tlist(i)
    fThresh = ThreshWav(f, T)


def exo2():
    """
    Test the effect of block thresholding on the image $f_0$ itself, for increasing value of $T$.
    (of course thresholding directly the image has no interest, this is just
    to vizualize the effect).
    """
    tlist = linspace(.3, .9, 4)
    for i in 1: length(tlist):
        T = tlist(i)
        imageplot(clamp(ThreshBlock(f, T)), ['T = ' num2str(T, 2)], 2, 2, i)


def exo3():
    """
    Display the evolution of the denoising SNR when $T$ varies.
    Store in |fBlock| the optimal denoising result.
    etrieve best
    """
    Tlist = linspace(.5, 2, 30)*sigma
    snr_stein = arrayfun(lambda T: snr(f0, ThreshWav(f, T)), Tlist)
    clf; plot(Tlist/ sigma, snr_stein(: ), 'LineWidth', 2)
    axis tight; set_label('T/ \sigma', 'SNR')
    [~, i] = max(snr_stein); T = Tlist(i)
    fBlock = ThreshWav(f, T)


def exo4():
    """
    Display the evolution of the denoising SNR when $T$ varies.
    Store in |fTI| the optimal denoising result.
    etrieve best
    """
    Tlist = linspace(.5, 2, 20)*sigma
    snr_stein = arrayfun(lambda T: snr(f0, ThreshWav(f, T)), Tlist)
    clf; plot(Tlist/ sigma, snr_stein(: ), 'LineWidth', 2)
    axis tight; set_label('T/ \sigma', 'SNR')
    [~, i] = max(snr_stein); T = Tlist(i)
    fTI = ThreshWav(f, T)


