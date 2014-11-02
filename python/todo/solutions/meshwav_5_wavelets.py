def exo1():
    """
    Plot the approximation curve error as a function of the number of
    coefficient.
    """


def exo2():
    """
    Perform denoising of spherical function by thresholding. Study the
    evolution of the optimal threshold as a function of the noise level.
    """


def exo3():
    """
    Display a dual
    wavelet that is used for the reconstruction
    by taking the inverse transform of a dirac.
    """
    nverts = size(vertex{end}, 2)
    i = 0
    for j  in  [J-3 J-2]:
        i = i + 1
        nj = size(vertex{j}, 2); nj1 = size(vertex{j + 1}, 2)
        sel = nj + 1: nj1-1
        d = sum(abs(vertex{end}(: , sel)))
        [tmp, k] = min(d); k = sel(k)
        fw2 = zeros(nverts, 1); fw2(k) = 1
        f2 = perform_wavelet_mesh_transform(vertex, face, fw2, -1, options)
        options.color = 'wavelets'
        options.use_color = 1
        options.rho = .6
        options.use_elevation = 1
        options.view_param = [104, -40]
        subplot(1, 2, i)
        plot_spherical_function(-vertex{end}, face{end}, f2, options); axis tight
        title(['Wavelet, scale ' num2str(j)])
    options.view_param = []


