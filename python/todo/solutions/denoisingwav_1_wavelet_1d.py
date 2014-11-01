def exo1():
    """
    Display the evolution of the denoising SNR
    when $T$ varies.
    Store in |fBest0| and |fBest1| the optimal denoising results.
    etrieve best
    """
    Tlist = linspace(.5, 6, 80)'*sigma
    snr0 = arrayfun(lambda T: snr(f0, Theta0W(f, T)), Tlist)
    snr1 = arrayfun(lambda T: snr(f0, Theta1W(f, T)), Tlist)
    clf; plot(Tlist/ sigma, [snr0 snr1], 'LineWidth', 2)
    axis tight; set_label('T/ \sigma', 'SNR')
    legend('Hard', 'Soft')
    [~, i] = max(snr0); fBest0 = Theta0W(f, Tlist(i))
    [~, i] = max(snr1); fBest1 = Theta1W(f, Tlist(i))


def exo2():
    """
    Display the evolution of the denoising SNR
    when $T$ varies.
    Store in |fBest0| and |fBest1| the optimal denoising results.
    etrieve best
    """
    Tlist = linspace(.5, 6, 80)'*sigma
    snr0 = arrayfun(lambda T: snr(f0, Theta0W(f, T)), Tlist)
    snr1 = arrayfun(lambda T: snr(f0, Theta1W(f, T)), Tlist)
    clf; plot(Tlist/ sigma, [snr0 snr1], 'LineWidth', 2)
    axis tight; set_label('T/ \sigma', 'SNR')
    legend('Hard TI', 'Soft TI')
    [~, i] = max(snr0); fBest0 = Theta0W(f, Tlist(i))
    [~, i] = max(snr1); fBest1 = Theta1W(f, Tlist(i))


