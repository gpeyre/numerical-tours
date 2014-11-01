def exo1():
    """
    Compute the true redundancy of the transform. Check that the transform
    is a tight frame (energy conservation).
    """
    P = 0
    for i in 1: length(S):
        P = P + length(S{i}(: ))
    disp(strcat(['True redundancy of the dictionary = ' num2str(P/ n) '.']))


def exo2():
    """
    Find the best threshold, that gives the smallest error.
    """
    Tlist = linspace(.4, 1.2, 15)*sigma
    err = []
    for i in 1: length(Tlist):
        ST = perform_thresholding(S, Tlist(i), 'soft')
        xT = perform_stft(ST, wlist, qlist, options)
        err(i) = snr(x0, xT)
    plot(Tlist/ sigma, err); axis('tight')
    set_label('T/ \sigma', 'SNR')


def exo3():
    """
    Perform the iterative thresholding by progressively decaying the value
    of |lambda| during the iterations, starting from |lambda=1.5*sigma| until
    |lambda=.5*sigma|. Retain the solution |xbp| together with the coefficients |Sbp|
    that provides the smallest
    error.
    """
    niter = 100
    lambda_max = 1.4*sigma
    lambda_min = .8*sigma
    lambda_list = linspace(lambda_max, lambda_min, niter)
    err = []
    for i in 1: niter:
        % progressbar(i, niter)
        lambda = lambda_list(i)
        % gradient
        r = x - x1
        Sr = perform_stft(r, wlist, qlist, options)
        S1 = cell_add(S1, Sr)
        % threshold
        S1 = perform_thresholding(S1, lambda, 'soft')
        % update the value of lambda to match noise
        x1 = perform_stft(S1, wlist, qlist, options)
        % lambda = lambda * sqrt(n)*sigma / norm(x-x1, 'fro')
        err(i) = snr(x0, x1)
        if i >1 & err(i) >max(err(1: i-1))
            xbp = x1; Sbp = S1
    plot(lambda_list/ sigma, err); axis('tight')
    set_label('\lambda/ \sigma', 'SNR')


def exo4():
    """
    Compare the source separation obtained by masking with a tight frame Gabor
    transform and with the coefficients computed by a basis pursuit
    sparsification process.
    """


