def exo1():
    """
    Recover an image from the transformed YUV representation $U_1$.
    """
    rgb2yuv = lambda f: applymat(f, T^(-1))
    f1 = rgb2yuv(U1)
    imageplot(f, 'Image', 1, 2, 1)
    imageplot(clamp(f1), 'Modified', 1, 2, 2)


def exo2():
    """
    Show the reconstructed image obtained by progressively rotating the hue component,
    i.e. replacing $H$ by $H + \theta$.
    his ugly code is for Scilab compatibility
    """
    a = lambda g: g(: , : , 2) .* cos(g(: , : , 1))
    b = lambda g: g(: , : , 2) .* sin(g(: , : , 1))
    c = lambda g: cat(3, g(: , : , 3), a(g), b(g))
    hsv12rgb = lambda g: applymat(c(g), T)
    theta = linspace(0, pi/ 2, 6)
    for i in 1: length(theta):
        g1 = g;  g1(: , : , 1) = g1(: , : , 1) + theta(i)
        imageplot(clamp(hsv12rgb(g1)), ['\theta = ' num2str(theta(i))], 2, 3, i)


def exo3():
    """
    Inverse the PCA transformation in order to retrive an image
    the modified PCA representation |U1|.
    isplay
    """
    pca1rgb = lambda g, V, m: applymat(g, V') + repmat(m, [n n 1])
    f1 = pca1rgb(g1, V, m)
    imageplot(f, 'f', 1, 2, 1)
    imageplot(clamp(f1), 'f_1', 1, 2, 2)


def exo4():
    """
    Compare (translation invariant) wavelet denoising of color image in the
    RGB and PCA color space (the PCA space should be estimated from the noisy image).
    """


def exo5():
    """
    For a complicated, non-homogenous color image, compute a local color
    space for each pixel, by performing a PCA over group of pixels. Use this
    adaptive color model to perform wavelet denoising.
    """


def exo6():
    """
    Compare wavelet domain image compression (quantization+coding) over the
    original RGB space and a more adapted global or local color space.
    """


