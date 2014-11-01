def exo1():
    """
    Perform the linear Fourier approximation with $M$ coefficients.
    Store the result in the variable |fM|.
    isplay
    """
    m = sqrt(M)
    F = fftshift(fft2(f))
    F1 = zeros(n0, n0)
    sel = (n0/ 2-m/ 2: n0/ 2 + m/ 2) + 1
    F1(sel, sel) = F(sel, sel)
    fM = real(ifft2(fftshift(F1)))
    imageplot(clamp(fM), ['Linear, Fourier, SNR = ' num2str(snr(f, fM), 4) 'dB'])


def exo2():
    """
    Compute the value of the threshold $T$ so that the number of
    coefficients is $M$. Display the corresponding approximation
    $f_M$.
    isplay
    """
    F = fft2(f)
    a = sort(abs(F(: )))
    if a(1) <a(n0*n0)
    	a = a(n0^2: -1: 1)
    T = a(M + 1)
    F = F .* (abs(F) >T)
    fM = real(ifft2(F))
    imageplot(clamp(fM), ['Non-linear, Fourier, SNR = ' num2str(snr(f, fM), 4) 'dB'])


def exo3():
    """
    Perform linear approximation with $M$ wavelet coefficients.
    isplay
    """
    q = sqrt(M)
    fw = perform_wavelet_transf(f, Jmin, + 1)
    fw1 = zeros(n0, n0)
    fw1(1: q, 1: q) = fw(1: q, 1: q)
    fM = perform_wavelet_transf(fw1, Jmin, -1)
    imageplot(clamp(fM), ['Linear, Wavelets, SNR = ' num2str(snr(f, fM), 4) 'dB'])


def exo4():
    """
    Perform non-linear approximation with $M$ wavelet coefficients by
    chosing the correct value for $T$. Store the result in the variable
    |fM|.
    isplay
    """
    fw = perform_wavelet_transf(f, Jmin, + 1)
    a = sort(abs(fw(: )))
    if a(1) <a(n0*n0)
    	a = a(n0^2: -1: 1)
    T = a(M + 1)
    fw1 = fw .* (abs(fw) >T)
    fM = perform_wavelet_transf(fw1, Jmin, -1)
    imageplot(clamp(fM), ['Non-linear, wavelets, SNR = ' num2str(snr(f, fM), 4) 'dB'])


