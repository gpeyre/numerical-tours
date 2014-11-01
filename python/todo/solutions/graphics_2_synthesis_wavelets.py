def exo1():
    """
    Iterate these two steps (spatial and wavelet histogram matching) until convergence to a stable step.
    pacial matching
    """
    niter = 4
    M1 = randn(n)
    for k in 1: niter:
        M1 = perform_hist_eq(M1, M)
        % wavelet matching
        MW1 = perform_wavelet_transf(M1, Jmin, + 1, options)
    for i in 1: size(MW, 3):
            MW1(: , : , i) = perform_hist_eq(MW1(: , : , i), MW(: , : , i))
        M1 = perform_wavelet_transf(MW1, Jmin, -1, options)
        % display
        imageplot(M1, strcat(['Iteration ' num2str(k)]), 2, 2, k)


def exo2():
    """
    Perform color texture synthesis with wavelets over the RGB space.
    isplay.
    """
    niter = 3
    M1 = randn(n, n, 3)
    for c in 1: 3:
        MW = perform_wavelet_transf(M(: , : , c), Jmin, + 1, options)
    for k in 1: niter:
            M1(: , : , c) = perform_hist_eq(M1(: , : , c), M(: , : , c))
            % wavelet matching
            MW1 = perform_wavelet_transf(M1(: , : , c), Jmin, + 1, options)
    for i in 1: size(MW, 3):
                MW1(: , : , i) = perform_hist_eq(MW1(: , : , i), MW(: , : , i))
            M1(: , : , c) = perform_wavelet_transf(MW1, Jmin, -1, options)
    imageplot(M, 'Exemplar', 1, 2, 1)
    imageplot(M1, 'Synthesized', 1, 2, 2)


def exo3():
    """
    Try with other color spaces, for instance PCA adapte space.
    """


def exo4():
    """
    Perform iteratively the randomized matching. Plot the decay of the
    mathing error.
    """
    M1 = randn(n, n, 3)
    niter1 = 200
    delta = []
    for k in 1: niter1:
        [U, R] = qr(randn(3))
        d = reshape(M, [n^2 3])*U
        d1 = reshape(M1, [n^2 3])*U
        delta(k) = 0
    for c in 1: 3:
            [tmp, I] = sort(d(: , c))
            [tmp, I1] = sort(d1(: , c))
            delta(k) = delta(k) + norm(d1(I1, c) - d(I, c))^2
            d1(I1, c) = d(I, c)
            % d1(: , c) = perform_hist_eq(d1(: , c), d(: , c))
        M1old = M1
        M1 = reshape(d1*U', [n n 3])
        % delta(k) = norm(M1(: )-M1old(: ))
    plot(log10(delta/ norm(M1(: ))), '.-')
    axis('tight')


def exo5():
    """
    Perform color texture synthesis with wavelets using this color histogram
    matching at each iteration.
    recompute the wavelet transform of the image
    isplay.
    """
    niter = 5
    niter1 = 5
    M1 = randn(n, n, 3)
    M1 = perform_color_matching(M1, M, niter1)
    for c in 1: 3:
        MW(: , : , : , c) = perform_wavelet_transf(M(: , : , c), Jmin, + 1, options)
    Msvg1 = {}; Msvg = {}
    for k in 1: niter:
        % wavelet matching
    for c in 1: 3:
            % transform
            MW1 = perform_wavelet_transf(M1(: , : , c), Jmin, + 1, options)
            % equalize
    for i in 1: size(MW1, 3):
                MW1(: , : , i) = perform_hist_eq(MW1(: , : , i), MW(: , : , i, c))
            M1(: , : , c) = perform_wavelet_transf(MW1, Jmin, -1, options)
        % spacial matching
        M1 = perform_color_matching(M1, M, niter1)
        Msvg{end + 1} = M1
    imageplot(M, 'Exemplar', 1, 2, 1)
    imageplot(M1, 'Synthesized', 1, 2, 2)


