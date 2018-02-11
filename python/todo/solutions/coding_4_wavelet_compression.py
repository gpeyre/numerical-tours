def exo1():
    """
    Compute the coefficients |MWT| obtained by thresholding at
    |T| the coefficients |MW|. Compute the coefficients |MWQ| obtained
    by quantizing with bin size |T| the same coefficients.
    Display them using the function |plot_wavelet|.
    hresholding approximation
    isplay
    """
    MWT = perform_thresholding(MW, T, 'hard')
    MWQ = perform_thresholding(MW, T, 'quantize')
    subplot(1, 2, 1)
    plot_wavelet(MWT, Jmin)
    title('Thresholded coefficients')
    subplot(1, 2, 2)
    plot_wavelet(MWT-MWQ, Jmin)
    title('Thresholded - Quantized')


def exo2():
    """
    Compare the effect of quantizing at |T=.2| and thresholding at |T=.2|
    the wavelet coefficients of an image.
    nverse transform
    rror
    isplay
    """
    MT = perform_wavelet_transf(MWT, Jmin, -1)
    MQ = perform_wavelet_transf(MWQ, Jmin, -1)
    eT = snr(M, MT)
    eQ = snr(M, MQ)
    imageplot(MT, strcat(['Thresholding, SNR = ' num2str(eT, 2)]), 1, 2, 1)
    imageplot(MT-MQ, strcat(['Thresholding - Approximating, SNR = + ' num2str(eT-eQ, 2)]), 1, 2, 2)


def exo3():
    """
    Compute a bin size |T0| to quantize the original |M| itself to obtained
    |MQ0| such that |norm(M-MQ,'fro')| is as close as possible to the error
    obtained with wavelet domain quantization.
    """
    Tlist = linspace(0.01, 0.3, 60)
    for i in 1: length(Tlist):
        err(i) = norm(M-perform_thresholding(M, Tlist(i), 'quantize'), 'fro')
    errW = norm(M-MQ, 'fro')
    [tmp, i] = min(abs(err-errW))
    T0 = Tlist(i)
    MQ0 = perform_thresholding(M, T0, 'quantize')
    disp(['Spatial quantization step T0 = ' num2str(T0, 2) '.'])
    imageplot(clamp(MQ), 'Wavelet quantized', 1, 2, 1)
    imageplot(clamp(MQ0), 'Spacial quantized', 1, 2, 2)


def exo4():
    """
    Compute the entropy lower bound for the quantized
    wavelet coefficients and for the quantized pixel values.
    Take care of |log2(0)| when |h(i)=0|.
    """
    h = h + 1e-10; h = h/ sum(h)
    E = -sum(h.*log2(h))
    h0 = h0 + 1e-10; h0 = h0/ sum(h0)
    E0 = -sum(h0.*log2(h0))
    disp(['Pixels entropy:  ' num2str(E0, 2)])
    disp(['Wavelet entropy: ' num2str(E, 2)])


def exo5():
    """
    Compute, for various threshold |T|, the number of bits per pixels |E(T)|
    of the quantized wavelet coefficients,
    and the wavelet decompression error |err(T)|, compute using SNR.
    Display the rate
    distortion curve |err| as a function of |E|.
    """
    Tlist = linspace(.03, .6, 20)
    err = []; nbits = []
    for i in 1: length(Tlist):
        T = Tlist(i)
        % quantize
        MWI = floor(abs(MW/ T)).*sign(MW)
        MWQ = sign(MWI) .* (abs(MWI) + .5) * T
        % inverse
        MQ = perform_wavelet_transf(MWQ, Jmin, -1)
        % error
        err(i) = snr(M, MQ)
        % bits
        nbits(i) = compute_entropy(MWI(: ))
    hh = plot(nbits, err); axis('tight')
    set_label('bpp', 'SNR')
    if using_matlab()
        set(hh, 'LineWidth', 2)


def exo6():
    """
    Extract the three fine scale wavelet coefficients (horizontal, vertical,
    diagonal directions) and quantize them, for instance with |T=.1|.
    Compute the entropy of the three sets together, and compute the entropy
    of each set.
    """
    Etot = compute_entropy([MWH MWV MWD])
    Ehor = compute_entropy(MWH)
    Ever = compute_entropy(MWV)
    Edia = compute_entropy(MWD)
    disp(['Entropy, all:  ' num2str(Etot, 3)])
    disp(['Entropy, hor:  ' num2str(Ehor, 3)])
    disp(['Entropy, vert: ' num2str(Ever, 3)])
    disp(['Entropy, diag: ' num2str(Edia, 3)])


def exo7():
    """
    Compare the number of bits needed to code all the wavelet coefficients
    together, and the number of bits needed to code independantly each scale
    of wavele coefficients for |Jmin=4<=j<=log2(n)-1| (and group together the
    remaining coefficients for |j<Jmin|).
    """
    Esep = 0
    Jmax = log2(n)-1; Jmin = 4
    for j  in  Jmax: -1: Jmin:
    for q in 1: 3:
            [selx, sely] = compute_quadsel(j, q)
            MWj = MWI(selx, sely)
            Esep = Esep + prod(size(MWj))*compute_entropy(MWj)
    Esep = Esep + prod(size(MWj))*compute_entropy(MWI(1: 2^j, 1: 2^j))
    Ewhole = compute_entropy(MWI)
    disp(['nb.bis, whole:    ' num2str(Ewhole, 3) ' bpp'])
    disp(['nb.bis, separate: ' num2str(Esep/ n^2, 3) ' bpp'])


def exo8():
    """
    Compute the rate distortion curve obtained by coding the coefficient
    separately through the scale, and compare with the rate distortion curve
    obtained by coding the coefficients as a whole.
    """
    nbits1 = []
    for i in 1: length(Tlist):
        T = Tlist(i)
        % quantize
        MWI = floor(abs(MW/ T)).*sign(MW)
        % bits    
        Esep = 0
        Jmax = log2(n)-1; Jmin = 4
    for j  in  Jmax: -1: Jmin:
    for q in 1: 3:
                [selx, sely] = compute_quadsel(j, q)
                MWj = MWI(selx, sely)
                Esep = Esep + prod(size(MWj))*compute_entropy(MWj)
        Esep = Esep + prod(size(MWj))*compute_entropy(MWI(1: 2^j, 1: 2^j))
        Ewhole = compute_entropy(MWI)
        nbits1(i) = Esep/ n^2
    hh = plot([nbits(: )'; nbits1(: )']', [err(: )'; err(: )']'); axis('tight')
    set_label('bpp', 'SNR')
    legend('Whole', 'Separate')
    if using_matlab()
        set(hh, 'LineWidth', 2)


