def exo1():
    """
    Compute the successive low pass approximations.
    """
    fj = f
    for j in 1: 4:
        fj = reshape(fj, [length(fj)/ 4 4])
        fj = mean(fj, 2)
        subplot(2, 2, j)
        options.face_vertex_color = fj
        plot_mesh(vertex{end-j}, face{end-j}, options)
        view(vv)
        colormap gray(256)
        lighting none


def exo2():
    """
    Compute the full wavelet transform, and check for orthogonality
    (conservation of energy).
    heck for orthogonality.
    """
    fw = f
    nj = length(f)
    for j in 1: J-1:
        fj = fw(1: nj)
        fj = reshape(fj, [nj/ 4 4])
        fj = fj*U
        fw(1: nj) = fj(: )
        nj = nj/ 4
    e = abs(norm(f)-norm(fw))/ norm(f)
    disp(['Orthogonality deviation (should be 0): ' num2str(e)])


def exo3():
    """
    Implement the backward spherical Haar transform (replace U by U' to perform the reconstruction),
    and check for perfect reconstruction.
    """
    nj = size(face{1}, 2)
    f1 = fw
    for j in 1: J-1:
        fj = f1(1: 4*nj)
        fj = reshape(fj, [nj 4])
        fj = fj*U'
        f1(1: 4*nj) = fj(: )
        nj = nj*4
    e = norm(f-f1)/ norm(f)
    disp(['Bijectivity deviation (should be 0): ' num2str(e)])


def exo4():
    """
    Perform Haar wavelet approximation with only 10% of the coefficients.
    orward transform
    hresholding
    ackward transform
    isplay
    """
    rho = .1
    rho = .05
    nj = length(f)
    fw = f
    for j in 1: J-1:
        fj = fw(1: nj)
        fj = reshape(fj, [nj/ 4 4])
        fj = fj*U
        fw(1: nj) = fj(: )
        nj = nj/ 4
    fwT = perform_thresholding(fw, round(rho*n), 'strict')
    nj = size(face{1}, 2)
    f1 = fwT
    for j in 1: J-1:
        fj = f1(1: 4*nj)
        fj = reshape(fj, [nj 4])
        fj = fj*U'
        f1(1: 4*nj) = fj(: )
        nj = nj*4
    options.face_vertex_color = clamp(f1)
    plot_mesh(vertex{end}, face{end}, options)
    view(vv)
    colormap gray(256)
    lighting none


def exo5():
    """
    Compare with the traditional 2D Haar approximation of |M|.
    """
    options.h = compute_wavelet_filter('Daubechies', 2)
    MW = perform_wavortho_transf(M, 0, + 1, options)
    MWT = perform_thresholding(MW(: ), round(rho*n), 'strict')
    MWT = reshape(MWT, size(M))
    M1 = perform_wavortho_transf(MWT, 0, -1, options)
    imageplot(clamp(crop(M1)))


def exo6():
    """
    Implement Spherical denoising using the Haar transform. Compare it with
    vertex-based lifting scheme denoising.
    """


