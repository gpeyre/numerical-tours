def exo1():
    """
    For a given $\lambda$, display the histogram of the repartition of
    the quadratic error $\norm{y-h(y)}^2$ and of $\text{SURE}(y)$.
    Compute these repartition using Monte-Carlo simulation (you need to
    generate lots of different realization of the noise $W$.
    Display in particular the location of the mean of these quantities.
    ean
    isplay
    """
    ntrials = 100; nlaunch = 20
    E0 = []; E = []
    for i in 1: nlaunch:
    	f = repmat(f0, [1 1 ntrials]) + sigma*randn(n, n, ntrials)
    	hf = h(f, lambda)
        %
        e = sum(sum((hf - repmat(f0, [1 1 ntrials])).^2))
        E0 = [E0; e(: )]
        %
        e = -N*sigma^2 + ...
            sum(sum((hf - f).^2)) + ...
            2 * sigma^2 * df(lambda); E = E(: )
        E = [E; e(: )]
    v_true = mean(E0)
    v_sure = mean(E)
    a = mean(E0) - 8*std(E0)
    b = mean(E0) + 8*std(E0)
    t = linspace(a, b, 31)
    mybar = lambda e: hist(e(e >a & e <b), t)
    subplot(2, 1, 1); hold on
    s = mybar(E0); bar(t, s); axis tight
    u = plot(v_true*[1 1], [0 max(s)], 'r'); set(u, 'LineWidth', 3)
    subplot(2, 1, 2); hold on
    s = mybar(E); bar(t, s); axis tight
    u = plot(v_sure*[1 1], [0 max(s)], 'r'); set(u, 'LineWidth', 3)


def exo2():
    """
    Compute, for a single realization $f=f_0+w$, the evolution
    of
    $$ E(\la) = \text{SURE}_\la(f) \qandq E_0(\lambda) = \norm{f-h_\la(f)}^2 $$
    as a function of $\lambda$.
    """
    f = f0 + sigma*randn(n)
    ntrials = 30
    lambda_list = linspace(.3, 2, ntrials)
    E0 = []; E = []
    for i in 1: ntrials:
        lambda = lambda_list(i)
        hf = h(f, lambda)
        E0(i) = norm(hf-f0, 'fro')^2
        E(i) = SURE(f, hf, lambda)
    plot(lambda_list, [E0; E])
    xlabel('\lambda'); legend('E_0', 'E')
    axis tight


def exo3():
    """
    Display the best denoising result $h_{\la^*}(f)$
    where
    $$\la^* = \uargmin{\la} \text{SURE}_\la(f) $$
    """
    [tmp, i] = min(E)
    lambda = lambda_list(i)
    hf = h(f, lambda)
    imageplot(clamp(hf), strcat(['SNR = ' num2str(snr(f0, hf), 3) 'dB']))


def exo4():
    """
    For a given $\lambda$, display the histogram of the repartition of
    the quadratic error $\norm{y-h(y)}^2$ and of $\text{SURE}(y)$.
    Compute these repartition using Monte-Carlo simulation (you need to
    generate lots of different realization of the noise $W$.
    Display in particular the location of the mean of these quantities.
    _Hint:_ you can do the computation directly over the wavelet domain,
    i.e. consider that the noise is added to the wavelet transform.
    ean
    isplay
    """
    f0w = W(f0)
    ntrials = 100; nlaunch = 20
    E0 = []; E = []
    dfw = lambda fw, lambda: sum(sum(abs(fw) >lambda))
    for i in 1: nlaunch:
        fw = repmat(f0w, [1 1 ntrials]) + sigma*randn(n, n, ntrials)
        hfw = S(fw, lambda)
        %
        e = sum(sum((hfw - repmat(f0w, [1 1 ntrials])).^2))
        E0 = [E0; e(: )]
        %
        e = -N*sigma^2 + ...
            sum(sum((hfw - fw).^2)) + ...
            2 * sigma^2 * dfw(fw, lambda); E = E(: )
        E = [E; e(: )]
    v_true = mean(E0)
    v_sure = mean(E)
    a = mean(E0) - 8*std(E0)
    b = mean(E0) + 8*std(E0)
    t = linspace(a, b, 31)
    mybar = lambda e: hist(e(e >a & e <b), t)
    subplot(2, 1, 1); hold on
    s = mybar(E0); bar(t, s); axis tight
    u = plot(v_true*[1 1], [0 max(s)], 'r'); set(u, 'LineWidth', 3)
    subplot(2, 1, 2); hold on
    s = mybar(E); bar(t, s); axis tight
    u = plot(v_sure*[1 1], [0 max(s)], 'r'); set(u, 'LineWidth', 3)


def exo5():
    """
    Compute, for a single realization $f=f_0+w$, the evolution
    of
    $$ E(\la) = \text{SURE}_\la(f) \qandq E_0(\lambda) = \norm{f-h_\la(f)}^2 $$
    as a function of $\lambda$.
    """
    f = f0 + sigma*randn(n)
    ntrials = 30
    lambda_list = linspace(.5, 2.5, ntrials)*sigma
    E0 = []; E = []
    for i in 1: ntrials:
        lambda = lambda_list(i)
        hf = h(f, lambda)
        E0(i) = norm(hf-f0, 'fro')^2
        E(i) = SURE(f, hf, lambda)
    plot(lambda_list/ sigma, [E0; E])
    xlabel('\lambda/ \sigma'); legend('E_0', 'E')
    axis tight


def exo6():
    """
    Display the best denoising result $h_{\la^*}(f)$
    where
    $$\la^* = \uargmin{\la} \text{SURE}_\la(f) $$
    """
    [tmp, i] = min(E)
    lambda = lambda_list(i)
    hf = h(f, lambda)
    imageplot(clamp(hf), strcat(['SNR = ' num2str(snr(f0, hf), 3) 'dB']))


