def exo1():
    """
    Compute the local Fourier transform around a point |t0| of |x|, which is the FFT (use the
    function |fft|) of the windowed signal |x.*h| where |h| is smooth
    windowing function located around |t0|. For instance you can use for |h|
    a Gaussian bump centered at |t0|. To center the FFT for display, use
    |fftshift|.
    enter for the Fourier analysis
    idth of the bump
    indow
    ft
    isplay signal
    isplay FFTs
    """
    t0 = n/ 4
    sigma = 128
    t = (1: n)'
    h = exp(-(t-t0).^2 ./ (2*sigma^2))
    xh = x.*h
    tau = 1e-3
    xh = xh(abs(xh) >tau)
    xf = fft(xh)
    subplot(2, 1, 1)
    plot(t, x); axis('tight')
    set_graphic_sizes([], 20)
    title('Signal x')
    subplot(2, 1, 2)
    plot(t, h); axis('tight')
    set_graphic_sizes([], 20)
    title('Window h')
    p = length(xf)
    subplot(2, 1, 1)
    plot(xh)
    set_graphic_sizes([], 20); axis('tight')
    title('Windowed signal (zoom)')
    subplot(2, 1, 2)
    plot(-p/ 2 + 1: p/ 2, fftshift(abs(xf)))
    set_graphic_sizes([], 20); axis('tight')
    title('Windowed FFT (zoom)')


def exo2():
    """
    A denoising is performed by hard or soft thresholding the STFT of the
    noisy signal. Compute the denosing SNR with both soft and hard
    thresholding, and compute the threshold that minimize the SNR. Remember that a soft thresholding
    should be approximately twice smaller than a hard thresholding. Check the
    result by listening. What can you conclude about the quality of the
    denoised signal ?
    etrieve best hard thresholding result
    isplay the error curves
    """
    Sn = perform_stft(xn, w, q, options)
    Tlist = linspace(.8, 2.5, 20)*sigma
    err_hard = []; err_soft = []
    for i in 1: length(Tlist):
        % soft thresholding
        SnT = perform_thresholding(Sn, Tlist(i)/ 2, 'soft')
        x1 = perform_stft(SnT, w, q, options)
        err_soft(i) = snr(x, x1)
        % hard thresholding
        SnT = perform_thresholding(Sn, Tlist(i), 'hard')
        x1 = perform_stft(SnT, w, q, options)
        err_hard(i) = snr(x, x1)
    [snr_hard, i] = max(err_hard)
    SnT = perform_thresholding(Sn, Tlist(i), 'hard')
    x1 = perform_stft(SnT, w, q, options)
    plot(Tlist/ sigma, [err_hard(: ) err_soft(: )])
    axis('tight')
    legend('Hard', 'Soft')
    set_graphic_sizes([], 20, 2)
    set_label('T/ sigma', 'SNR')


def exo3():
    """
    Display and hear the results. What do you notice ?
    isplay
    ear
    """
    subplot(2, 1, 1)
    plot(xn); axis([1 n -1.2 1.2])
    set_graphic_sizes([], 20)
    title(strcat(['Noisy signal, SNR = ', num2str(snr(x, xn), 4), 'dB']))
    subplot(2, 1, 2)
    plot(x1); axis([1 n -1.2 1.2])
    set_graphic_sizes([], 20)
    title(strcat(['Denoised signal, SNR = ', num2str(snr_hard, 4), 'dB']))
    sound(x1, fs)


def exo4():
    """
    Trie for various block sizes and report the best results.
     progressbar(k,length(bsX(:)));
    """
    Sn = perform_stft(xn, w, q, options)
    Tlist = linspace(1, 3, 6)/ 2*sigma
    snr_block = []
    bslist = 1: 4
    [bsY, bsX] = meshgrid(bslist, bslist)
    for k in 1: length(bsX(: )):
        options.block_size = [bsX(k) bsY(k)]
        err = []
    for i in 1: length(Tlist):
            % hard thresholding
            SnT = perform_thresholding(Sn, Tlist(i), 'block', options)
            x1 = perform_stft(SnT, w, q, options)
            err(i) = snr(x, x1)
        % retrieve best hard thresholding result
        [snr_block(k), t] = max(err)
        if t = =1 | t = =length(Tlist)
            warning('Out of bound reached')
        Topt(k) = Tlist(t)
    snr_block = reshape(snr_block, size(bsX))
    imageplot(snr_block, 'SNR for several block sizes')
    set_label('X block size', 'Y block size')


