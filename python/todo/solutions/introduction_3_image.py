def exo1():
    """
    To avoid boundary artifacts and estimate really the frequency content of
    the image (and not of the artifacts!), one needs to multiply |M| by a
    smooth windowing function |h| and compute |fft2(M.*h)|. Use a sine
    windowing function. Can you interpret the resulting filter ?
    ompute kernel h
    ompute FFT
    isplay
    """
    t = linspace(-pi(), pi(), n)
    h = (cos(t) + 1)/ 2
    h = h'*h
    Mf = fft2(M.*h)
    Lf = fftshift(log(abs(Mf) + 1e-1))
    imageplot(M.*h, 'Image', 1, 2, 1)
    imageplot(Lf, 'Fourier transform', 1, 2, 2)


def exo2():
    """
    Perform low pass filtering by removing the high frequencies of the
    spectrum. What do you oberve ?
    isplay
    """
    k = round(.8*n); k = round(k/ 2)*2; % even number
    Mf = fft2(M)
    Mf(n/ 2-k/ 2 + 2: n/ 2 + k/ 2, n/ 2-k/ 2 + 2: n/ 2 + k/ 2) = 0
    Mh = real(ifft2(Mf))
    imageplot(crop(M), 'Image', 1, 2, 1)
    imageplot(clamp(crop(Mh)), 'Low pass filtered', 1, 2, 2)


