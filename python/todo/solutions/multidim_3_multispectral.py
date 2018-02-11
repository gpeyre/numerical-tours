def exo1():
    """
    Implement the inverse transform to recover an approximation |M1| from
    the coefficients |UWT|.
    proximation error
    isplay
    """
    M1 = UWT
    for i in 1: p:
        M1(: , : , i) = perform_wavelet_transf(M1(: , : , i), Jmin, -1)
    M1 = reshape(M1, [n*n p])'
    M1 = idct(M1)
    M1 = reshape(M1', [n n p])
    e = snr(M, M1)
    imageplot(M(: , : , rgbsel), 'Original', 1, 2, 1)
    imageplot(clamp(M1(: , : , rgbsel)), strcat(['Approximated, SNR = ' num2str(e, 2) 'dB']), 1, 2, 2)


def exo2():
    """
    Compare the approximation error (both in term of SNR and visually)
    of a multispectral image with a 3D Haar basis and with a tensor product
    of a 2D Haar and a DCT.
    """


def exo3():
    """
    Compare the denoising (both in term of SNR and visually) of a multispectral image with an independant
    thresholding of each channel within a translation invariant 2D wavelet basis,
    and with a thresholding of the DCT/invariant wavelet representation. For
    each method, compute the optimal threshold value.
    """


