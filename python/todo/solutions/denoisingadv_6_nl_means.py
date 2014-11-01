def exo1():
    """
    Compute the denoising result for several values of $\tau$ in order to
    determine the optimal denoising that minimizes $\norm{\tilde f - f_0}$.
    """
    ntests = 4
    tau_list = linspace(.003, .025, ntests)
    e0 = -Inf
    for k in 1: ntests:
        tau = tau_list(k)
        f1 = NLmeans(tau)
        e = snr(f0, f1)
        if e >e0
            fNL = f1
            e0 = e
        imageplot(clamp(f1), ['SNR = ' num2str(e, 4) 'dB'], 2, 2, k)


def exo2():
    """
    Explore the influence of the $q$ and $w$ parameters.
    """


