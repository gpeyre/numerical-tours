def exo1():
    """
    Display a filtering $\Gg_\si(f_0)$ with increasing size $\si$.
    """
    slist = [.5 3 6 10]
    for i in 1: 4:
        s = slist(i)
        imageplot(Filter(f0, s), ['\sigma = ' num2str(s)], 2, 2, i)


def exo2():
    """
    Display several weights $ W_{v_i} $.
    """
    vlist = round([.1 .3 .6 .9]*p)
    for i in 1: 4:
        v = vlist(i)
        imageplot(W(: , : , v), ['v_i = ' num2str((v-1)/ (p-1), 2)], 2, 2, i)


def exo3():
    """
    Display several stacks $F_{v_i}$.
    """
    vlist = round([.1 .3 .6 .9]*p)
    for i in 1: 4:
        v = vlist(i)
        imageplot(F(: , : , v), ['v_i = ' num2str((v-1)/ (p-1), 2)], 2, 2, i)


def exo4():
    """
    Compare nearest-neighbor and linear destacking.
    """
    fNN = bilateral_nn(f0, sx, sv)
    fLin = bilateral_lin(f0, sx, sv)
    c = [.5 .4]*n; q = 200
    imageplot(crop(fNN, q, c), 'Nearest neighbor', 1, 2, 1)
    imageplot(crop(fLin, q, c), 'Linear', 1, 2, 2)


def exo5():
    """
    Study the influence of $\sigma_x$ on the filter, for a fixed
    $\sigma_v=0.2$.
    """
    sv = .2
    slist = linspace(1, 8*2, 4)
    for i in 1: 4:
        sx = slist(i)
        imageplot(bilateral_lin(f0, sx, sv), ['\sigma_x = ' num2str(sx, 2)], 2, 2, i)


def exo6():
    """
    Study the influence of $\sigma_v$ on the filter, for a fixed
    $\sigma_x=8$.
    """
    sx = 4*2
    slist = linspace(.05, .4, 4)
    for i in 1: 4:
        sv = slist(i)
        imageplot(bilateral_lin(f0, sx, sv), ['\sigma_v = ' num2str(sv, 2)], 2, 2, i)


def exo7():
    """
    Compute the optimal parameter $(\sigma_x,\sigma_v)$ to maximize the
    SNR between $f_0$ and the filtered image. Record the optimal denoising
    result in |fOpt|.
    """
    if 0
    qx = 10; qv = 10
    sxlist = linspace(.5, 2*1.5, qx)
    svlist = linspace(.05, .3, qv)
    err = []
    for ix in 1: qx:
    for iv in 1: qv:
            sx = sxlist(ix); sv = svlist(iv)
            f1 = bilateral_lin(f, sx, sv)
            err(ix, iv) = snr(f0, f1)
            if err(ix, iv) = =max(err(: ))
                fOpt = f1
    surf(sxlist, svlist, err')
    axis('tight')
    view(60, 40)
    colormap(jet(256))
    xlabel('\sigma_x')
    ylabel('\sigma_v')
    else
        fOpt = f


def exo8():
    """
    Compare with translation invariant wavelet thresholding.
    """
    Jmin = 3
    options.ti = 1
    fW = perform_wavelet_transf(f, Jmin, + 1, options)
    T = 2.8*mu
    fW = fW .* (abs(fW) >T)
    fWav = perform_wavelet_transf(fW, Jmin, -1, options)
    imageplot(clamp(fWav), ['Wavelets, SNR = ' num2str(snr(f0, fWav), 3) 'dB'])


def exo9():
    """
    Perform detail boosting by enhancing the detail layer.
    For instance use various non-linear remapping of the intensities such as
    $ f_1 + \ga r^\al $ for some values of $\ga \geq 1$ and $\al>0$.
    """
    glist = [1 1.2 2 4]
    for i in 1: 4:
        gamma = glist(i)
        fr = f1 + gamma*r
        imageplot(clamp(fr), ['\gamma = ' num2str(gamma)], 2, 2, i)


def exo10():
    """
    Extend the bilateral filter for color images.
    """


def exo11():
    """
    Try several tone mapping operators, such as for instance
    $$ \phi_1(t)=\frac{t}{t+\epsilon} \qandq \phi_2(t) = \log(t+\epsilon) $$
    """
    imageplot(color_recompose(log(fV + 1e-5)) , '\phi_1', 1, 2, 1)
    imageplot(color_recompose(fV./ (fV + 1e-4)), '\phi_2', 1, 2, 2)


def exo12():
    """
    Compute the tone mapped image using $\tilde f_V$.
    Test with several value of $\ga,\epsilon, \si_x,\si_v$.
    """
    remap = lambda FV, FV1, gamma: exp((gamma*FV1 + FV-FV1) * (b-a) + a) - epsilon
    gamma = .1
    FVmapped = remap(FV, FV1, gamma)
    imageplot(color_recompose(FVmapped), ['\gamma = ' num2str(gamma)])


