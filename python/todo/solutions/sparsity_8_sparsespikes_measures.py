def exo1():
    """
    You can see that the dual certificate $\abs{\eta_\la}$ saturate to
    $+1$ or $-1$ in region that are far away from the initial support
    $x_0$. This means either that the noise is too large for the method to
    successfully estimate robustly this support, or that $\la$ was not
    chosen large enough. Explore different value of noise level $\norm{w}$ and $\la$
    to determine empirically the signal/noise/$\la$ regime where the support is
    sucessfully estimated.
    """


def exo2():
    """
    Study the evolution of the pre-certificate as the separation between the spikes diminishes.
    When is it the case that $\eta_0=\eta_V$? When is it the case that $\mu_0$ is the solution of $\Pp_0(\Phi \mu_0)$ ?
    """
    dList = [.4 .6 .8 1.2]
    for i in 1: length(dList):
        delta = dList(i)/ fc
        x0 = [.5-delta .5 .5 + delta]'
        %
        w = ones(N, 1)
        Gamma = []
    for k in 0: d:
            Gamma = [Gamma, diag(w) * Fourier(fc, x0)]
            % derivate the filter
            w = w .* 2i*pi .* (-fc: fc)'
        %
        pV = pinv(Gamma') * [sign(a0); zeros(n, 1)]
        etaV = PhiS(fc, u, pV)
        %
        subplot(2, 2, i)
        hold on
        stem(x0, sign(a0), 'k.--', 'MarkerSize', ms, 'LineWidth', lw)
        plot([0 1],  [1 1], 'k--', 'LineWidth', lw)
        plot([0 1], -[1 1], 'k--', 'LineWidth', lw)
        plot(u, etaV, 'b', 'LineWidth', lw)
        axis([0 1 -1.4 1.4])
        set(gca, 'XTick', [], 'YTick', [-1 1])
        box on
        title(['\delta = ' num2str(dList(i)) '/ f_c'])


