def exo1():
    """
    Compute the best threshold to minimize the denoising error in curvelets.
    Call |Mcurv| the optimal denoising.
    """
    Tlist = linspace(.8, 1.2, 15)*sigma
    MW = perform_curvelet_transform(M, options)
    err = []
    for i in 1: length(Tlist):
        MWT = perform_thresholding(MW, Tlist(i), 'hard')
        M1 = perform_curvelet_transform(MWT, options)
        err(end + 1) = snr(M0, M1)
    h = plot(Tlist/ sigma, err)
    set(h, 'LineWidth', 2)
    axis('tight')
    [tmp, i] = max(err)
    MWT = perform_thresholding(MW, Tlist(i), 'hard')
    Mcurv = perform_curvelet_transform(MWT, options)


def exo2():
    """
    Perform cycle spinning to enhance the recovery error.
    """
    T = Tlist(i)
    s = 4
    [dY, dX] = meshgrid(0: s-1, 0: s-1)
    Mcurv = zeros(n)
    for i in 1: s^2:
        Ms = circshift(M, [dX(i) dY(i)])
        MW = perform_curvelet_transform(Ms, options)
        MWT = perform_thresholding(MW, T, 'hard')
        Ms = perform_curvelet_transform(MWT, options)
        Ms = circshift(Ms, -[dX(i) dY(i)])
        Mcurv = (1-1/ i)*Mcurv + 1/ i*Ms


def exo3():
    """
    Compare with translation invariant hard thresholding.
    """
    options.ti = 1
    Jmin = 3
    T = 2.8*sigma
    MW = perform_wavelet_transf(M, Jmin, + 1, options)
    MWT = perform_thresholding(MW, T, 'hard')
    Mwav = perform_wavelet_transf(MWT, Jmin, -1, options)
    imageplot(clamp(Mwav), ['Wavelets, SNR = ' num2str(snr(M0, Mwav), 3) 'dB'], 1, 2, 1)
    imageplot(clamp(Mcurv), ['Curvelet, SNR = ' num2str(snr(M0, Mcurv), 3) 'dB'], 1, 2, 2)


def exo4():
    """
    Applies curvelet iterative thresholding to solve an inverse problem such
    as inpainting, deconvolution or compressed sending.
    """


