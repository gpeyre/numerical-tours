def exo1():
    """
    Try for various Gaussian variance |mu| to compute the denoising |xh|.
    Compute, in an oracle manner, the best variance |muopt| by computing the
    residual error |snr(x0,xh)|.
    isplay blurring for various mu
    ompute the error for many mu
    etrieve the best denoising result
    """
    mulist = linspace(.1, 3.5, 31)
    err =  []
    for i in 1: length(mulist):
        mu = mulist(i)
        % compute the filter   
        h = exp(-(t.^2)/ (2*mu^2))
        h = h/ sum(h(: ))
        % perform blurring
        xh = real(ifft(fft(x) .* fft(fftshift(h))))
        err(i) = snr(x0, xh)
    plot(mulist, err, '.-'); axis('tight')
    set_label('\mu', 'SNR')
    [snr_opt, I] = max(err)
    muopt = mulist(I)
    disp(strcat(['The optimal smoothing width is ' num2str(muopt) ' pixels, SNR = ' num2str(snr_opt) 'dB.']))


def exo2():
    """
    Try for various Gaussian variance to compute the denoising |Mh|.
    Compute, in an oracle manner, the best variance |muopt| by computing the
    residual error |snr(M0,Mh)|.
    ow compute the error for many mu
    etrieve the best denoising result
    """
    mulist = linspace(.3, 6, 31)
    err =  []
    for i in 1: length(mulist):
        mu = mulist(i)
        Mh = perform_blurring(M, mu, options)
        err(i) = snr(M0, Mh)
    plot(mulist, err, '.-'); axis('tight')
    set_label('\mu', 'SNR')
    [snr_opt, I] = max(err)
    muopt = mulist(I)
    disp(strcat(['The optimal smoothing width is ' num2str(muopt) ' pixels, SNR = ' num2str(snr_opt) 'dB.']))


