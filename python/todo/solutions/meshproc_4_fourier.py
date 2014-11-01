def exo1():
    """
    Show the smoothed mesh for an increasing number of Fourier frequencies |nb|.
    """
    nblist = round(linspace(10, nb, 4))
    for i in 1: length(nblist):
        V = U(: , 1: nblist(i))
        subplot(2, 2, i)
        plot_mesh((vertex*V)*V', faces)


def exo2():
    """
    Compute a best |m|-term non-linear approximation whith |m=.1*n|, by
    hard thresholding the Fourier coefficients using the correct threshold.
    Compare with linear |m| term approximation (use |m/3| coefficient for each
    coordinate X/Y/Z).
    on linear
    inear
    isplay
    """
    m = round(.1*n/ 3)*3
    pvertex = vertex*U
    pvertexN = perform_thresholding(pvertex, m, 'largest')
    vertexN = pvertexN*U'
    pvertexL = pvertex
    pvertexL(: , m/ 3 + 1: n) = 0
    vertexL = pvertexL*U'
    subplot(1, 2, 1)
    plot_mesh(vertexL, faces)
    subplot(1, 2, 2)
    plot_mesh(vertexN, faces)
    disp(['Linear:     SNR = ' num2str(snr(vertex, vertexL), 3) 'dB'])
    disp(['Non-linear: SNR = ' num2str(snr(vertex, vertexN), 3) 'dB'])


def exo3():
    """
    Compare the rate-distortion curve (log of error as a function of the
    log of the number of coefficients) for linear and non-linear approximation.
    onlinear
    inear
    ormalize
    isplay
    """
    v = sort(abs(pvertex(: )).^2, 'ascend')
    errN = reverse(cumsum(v))
    v = reverse(abs(pvertex(: )).^2)
    errL = reverse(cumsum(v))
    errL = errL/ errL(1)
    errN = errN/ errN(1)
    close; clf
    plot(log10(1: 3*n), log10([errL errN])); axis('tight')
    axis([0 log10(2*n) -5 0])
    legend('Linear', 'Non-linear')


def exo4():
    """
    Perform the compression for several quantization steps |T|
    and display the rate distortion curve showing the SNR
    as a function of the number of bits.
    lot
    """
    ntests = 50
    Tlist = linspace(.1, 3, ntests)
    err = []
    for i in 1: ntests:
        T = Tlist(i)
        % decoding
        pvertexI = floor(abs(pvertex/ T)).*sign(pvertex)
        pvertexQ = sign(pvertexI) .* (abs(pvertexI) + .5) * T
        vertex1 = pvertexQ*U'
        % entropic
        t = min(pvertexI(: )): max(pvertexI(: ))
        h = hist(pvertexI(: ), t)
        h = max(h, 1e-10); h = h/ sum(h)
        E = -sum(log2(h).*h)
        % recode
        nbits(i) = 3*E
        err(i) = snr(vertex, vertex1)
    plot(nbits, err); axis('tight')
    set_label('nb.bits', 'SNR')


