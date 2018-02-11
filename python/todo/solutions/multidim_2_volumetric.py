def exo1():
    """
    Implement the forward wavelet transform by iteratively applying these
    transform steps to the low pass residual.
    nitialize the transform
    """
    MW = M
    for j in 1: log2(n):
        p = n/ 2^(j-1)
        sel = 1: p
        % average/ difference along X
        MW(sel, sel, sel) = cat3(1, (MW(1: 2: p, sel, sel) + MW(2: 2: p, sel, sel))/ sqrt(2), (MW(1: 2: p, sel, sel)-MW(2: 2: p, sel, sel))/ sqrt(2))
        % average/ difference along Y
        MW(sel, sel, sel) = cat3(2, (MW(sel, 1: 2: p, sel) + MW(sel, 2: 2: p, sel))/ sqrt(2), (MW(sel, 1: 2: p, sel)-MW(sel, 2: 2: p, sel))/ sqrt(2))
        % average/ difference along Z
        MW(sel, sel, sel) = cat3(3, (MW(sel, sel, 1: 2: p) + MW(sel, sel, 2: 2: p))/ sqrt(2), (MW(sel, sel, 1: 2: p)-MW(sel, sel, 2: 2: p))/ sqrt(2))


def exo2():
    """
    Implement the backward transform to compute an approximation |M1| from
    the coefficients |MWT|.
    """
    M1 = MWT
    for j in log2(n): -1: 1:
        p = n/ 2^(j)
        sel = 1: p
        sel1 = 1: 2*p
        selw = p + 1: 2*p
        % average/ difference along X
        A = M1(sel,  sel1, sel1)
        D = M1(selw, sel1, sel1)
        M1(1: 2: 2*p, sel1, sel1) = (A + D)/ sqrt(2)
        M1(2: 2: 2*p, sel1, sel1) = (A-D)/ sqrt(2)
        % average/ difference along Y
        A = M1(sel1, sel,  sel1)
        D = M1(sel1, selw, sel1)
        M1(sel1, 1: 2: 2*p, sel1) = (A + D)/ sqrt(2)
        M1(sel1, 2: 2: 2*p, sel1) = (A-D)/ sqrt(2)
        % average/ difference along Z
        A = M1(sel1, sel1, sel)
        D = M1(sel1, sel1, selw)
        M1(sel1, sel1, 1: 2: 2*p) = (A + D)/ sqrt(2)
        M1(sel1, sel1, 2: 2: 2*p) = (A-D)/ sqrt(2)


def exo3():
    """
    Select the optimal blurring width |s| to reach the smallest possible
    SNR. Keep the optimal denoising |Mblur|
    """
    ntests = 20
    slist = linspace(.01, 1.5, ntests)
    err = []
    for i in 1: ntests:
        h = exp(-(X.^2 + Y.^2 + Z.^2)/ (2*slist(i)^2))
        h = h/ sum(h(: ))
        Mh = real(ifftn(fftn(Mnoisy) .* fftn(fftshift(h))))
        err(i) = snr(M, Mh)
        if i >1 && err(i) >max(err(1: i-1))
            Mblur = Mh
    plot(slist, err, '.-')
    axis('tight')
    set_label('s', 'SNR')


def exo4():
    """
    Perforn Wavelet denoising by thresholding the wavelet coefficients of
    Mnoisy. Test both hard thresholding and soft thresholding to determine
    the optimal threshold and the corresponding SNR.
    Record the optimal result |Mwav|.
    """
    MW = perform_haar_transf(Mnoisy, 1, + 1)
    Tlist = linspace(1, 4, 20)*sigma
    err_hard = []; err_soft = []
    for i in 1: length(Tlist):
        MWT = perform_thresholding(MW, Tlist(i), 'hard')
        M1 = perform_haar_transf(MWT, 1, -1)
        err_hard(i) = snr(M, M1)
        MWT = perform_thresholding(MW, Tlist(i), 'soft')
        M1 = perform_haar_transf(MWT, 1, -1)
        err_soft(i) = snr(M, M1)
        if i >1 & err_soft(i) >max(err_soft(1: i-1))
            Mwav = M1
    plot(Tlist/ sigma, [err_hard; err_soft]', '.-')
    axis('tight')
    set_label('T/ sigma', 'SNR')
    legend('hard', 'soft')


def exo5():
    """
    Implement cycle spinning hard thresholding with |T=3*sigma|.
    """
    T = 3*sigma
    w = 4
    [dX, dY, dZ] = ndgrid(0: w-1, 0: w-1, 0: w-1)
    Mspin = zeros(n, n, n)
    for i in 1: w^3:
        MnoisyC = circshift(Mnoisy, [dX(i) dY(i) dZ(i)])
        % denoise
        MW = perform_haar_transf(MnoisyC, 1, + 1)
        MWT = perform_thresholding(MW, T, 'hard')
        M1 = perform_haar_transf(MWT, 1, -1)
        % back
        M1 = circshift(M1, -[dX(i) dY(i) dZ(i)])
        Mspin = Mspin*(i-1)/ i + M1/ i


def exo6():
    """
    Implement the full 3D forward wavelet transform by applying these steps
    for decaying scales |j| toward 0.
    """
    Jmin = 0
    options.h = h
    MW = perform_wavortho_transf(M, Jmin, + 1, options)


def exo7():
    """
    Implement the full 3D backward wavelet transform by applying these steps
    for increasing scales |j|.
    """
    M1 = perform_wavortho_transf(MWT, Jmin, -1, options)


def exo8():
    """
    Implement denoising by soft and hard thresholding Daubechies wavelet
    coefficients.
    """
    MW = perform_wavortho_transf(Mnoisy, 1, + 1, options)
    Tlist = linspace(1, 4, 10)*sigma
    err_hard = []; err_soft = []
    for i in 1: length(Tlist):
        MWT = perform_thresholding(MW, Tlist(i), 'hard')
        M1 = perform_wavortho_transf(MWT, 1, -1, options)
        err_hard(i) = snr(M, M1)
        MWT = perform_thresholding(MW, Tlist(i), 'soft')
        M1 = perform_haar_transf(MWT, 1, -1, options)
        err_soft(i) = snr(M, M1)
        if i >1 & err_soft(i) >max(err_soft(1: i-1))
            Mwav = M1
    plot(Tlist/ sigma, [err_hard; err_soft]', '.-')
    axis('tight')
    set_label('T/ sigma', 'SNR')
    legend('hard', 'soft')


def exo9():
    """
    Implement cycle spinning hard thresholding with Daubechies wavelets with |T=3*sigma|.
    """
    T = 3*sigma
    w = 4
    [dX, dY, dZ] = ndgrid(0: w-1, 0: w-1, 0: w-1)
    Mspin = zeros(n, n, n)
    for i in 1: w^3:
        MnoisyC = circshift(Mnoisy, [dX(i) dY(i) dZ(i)])
        % denoise
        MW = perform_wavortho_transf(MnoisyC, 1, + 1, options)
        MWT = perform_thresholding(MW, T, 'hard')
        M1 = perform_wavortho_transf(MWT, 1, -1, options)
        % back
        M1 = circshift(M1, -[dX(i) dY(i) dZ(i)])
        Mspin = Mspin*(i-1)/ i + M1/ i


