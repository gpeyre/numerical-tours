def exo1():
    """
    Compute a best $M$-term approximation in the Fourier basis of $f$, for
    $M \in \{N/100, N/20\}$. Compute the approximation using
    a well chosen hard threshold value $T$.
    """
    Mlist = round([.01 .05]*n^2)
    for i in 1: length(Mlist):
        M = Mlist(i)
        fFT = perform_thresholding(fF, M, 'largest')
        fM = real(ifft2(fFT)*n)
        imageplot(clamp(fM), ['M/ N = ' num2str(M/ n^2, 2) ', SNR = ' num2str(snr(f, fM), 3) 'dB'], 1, 2, i)


def exo2():
    """
    Compute and display in log scales the ordered coefficients $c_R$.
    _Hint:_ a discrete primitive can be computed using the function
    |cumsum|.
    """
    cR = sort(abs(fF(: )))
    if cR(n^2) >cR(1)
        cR = reverse(cR); % be sure it is in reverse order
    lw = 2
    h = plot(log10(cR))
    if using_matlab()  
        set(h, 'LineWidth', lw)
    axis('tight')


def exo3():
    """
    Compute and display in log-scale the non-linear approximation
    error $\epsilon[M]^2$.
    Store the values of $\epsilon[M]^2$ in a vector |err_fft|.
    """
    err_fft = max(norm(f(: ))^2 - cumsum(cR.^2), 1e-10)
    h = plot(log10(err_fft / norm(f(: ))^2))
    if using_matlab()  
        set(h, 'LineWidth', lw)
    ax = [1 n^2/ 50 -2.35 0]
    axis(ax)
    title('log_{10}(\epsilon^2[M]/ ||f||^2)')


def exo4():
    """
    Compute a best $M$-term approximation in the wavelet basis of $f$, for
    $M \in \{N/100, N/20\}$. Compute the approximation using
    a well chosen hard threshold value $T$.
    Note that the inverse wavelet transform is obtained by replacing the +1
    by a -1 in the definition of the transform.
    """
    for i in 1: length(Mlist):
        M = Mlist(i)
        fWT = perform_thresholding(fW, M, 'largest')
        fM = perform_wavortho_transf(fWT, Jmin, -1, options)
        imageplot(clamp(fM), ['M/ N = ' num2str(M/ n^2, 2) ', SNR = ' num2str(snr(f, fM), 3) 'dB'], 1, 2, i)


def exo5():
    """
    Compute and display in log-scale the non-linear approximation
    error $\epsilon[M]^2$.
    Compares the Fourier and wavelets approximations.
    Store the values of $\epsilon[M]^2$ in a vector |err_wav|.
    """
    cR = sort(abs(fW(: )))
    if cR(n^2) >cR(1)
        cR = reverse(cR); % be sure it is in reverse order
    err_wav = max(norm(f(: ))^2 - cumsum(cR.^2), 1e-10)
    clf; hold('on')
    h = plot(log10(err_fft / norm(f(: ))^2), 'r')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    h = plot(log10(err_wav / norm(f(: ))^2), 'b')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    axis(ax)
    title('log_{10}(\epsilon^2[M]/ ||f||^2)')
    legend('Fourier', 'Wavelets')
    if using_matlab()  
        box('on')


def exo6():
    """
    Compute a best $M$-term approximation in the wavelet basis of $f$, for
    $M \in \{N/100, N/20\}$. Compute the approximation using
    a well chosen hard threshold value $T$. Note that the inverse DCT
    transform is obtained with the function |idct2|.
    """
    for i in 1: length(Mlist):
        M = Mlist(i)
        fCT = perform_thresholding(fC, M, 'largest')
        fM = idct2(fCT)
        imageplot(clamp(fM), ['M/ N = ' num2str(M/ n^2, 2) ', SNR = ' num2str(snr(f, fM), 3) 'dB'], 1, 2, i)


def exo7():
    """
    Compute and display in log-scale the non-linear approximation
    error $\epsilon[M]^2$.
    Compares the Fourier and DCT approximations.
    Store the values of $\epsilon[M]^2$ in a vector |err_dct|.
    """
    cR = sort(abs(fC(: )))
    if cR(n^2) >cR(1)
        cR = reverse(cR); % be sure it is in reverse order
    err_dct = max(norm(f(: ))^2 - cumsum(cR.^2), 1e-10)
    clf; hold('on')
    h = plot(log10(err_fft / norm(f(: ))^2), 'r')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    h = plot(log10(err_dct / norm(f(: ))^2), 'b')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    axis(ax)
    title('log_{10}(\epsilon^2[M]/ ||f||^2)')
    if using_matlab()  
        box('on')
    legend('Fourier', 'DCT')


def exo8():
    """
    Compute the local DCT transform |fL| by transforming each patch.
    """
    fL = zeros(n, n)
    for i in 1: n/ w:
    for j in 1: n/ w:
            seli = (i-1)*w + 1: i*w
            selj = (j-1)*w + 1: j*w
            fL(seli, selj) = dct2(f(seli, selj))


def exo9():
    """
    Compute the inverse local DCT transform  of the coefficients |fL| by inverse
    transforming each patch using  the function |idct2|.
    """
    f1 = fL
    for i in 1: n/ w:
    for j in 1: n/ w:
            seli = (i-1)*w + 1: i*w
            selj = (j-1)*w + 1: j*w
            f1(seli, selj) = idct2(f1(seli, selj))
    disp(strcat((['Error |f-f1|/ |f| = ' num2str(norm(f(: )-f1(: ))/ norm(f(: )))])))


def exo10():
    """
    Compute a few best |m|-term approximations in the Local DCT basis of
    |f|.
    """
    for u in 1: length(Mlist):
        M = Mlist(u)
        fLT = perform_thresholding(fL, M, 'largest')
        % inverse
        fM = fLT
    for i in 1: n/ w:
    for j in 1: n/ w:
                seli = (i-1)*w + 1: i*w
                selj = (j-1)*w + 1: j*w
                fM(seli, selj) = idct2(fM(seli, selj))
        % display
        imageplot(clamp(fM), ['M/ n^2 = ' num2str(M/ n^2, 2) ', SNR = ' num2str(snr(f, fM), 3) 'dB'], 1, 2, u)


def exo11():
    """
    Compute and display in log-scale the non-linear approximation
    error $\epsilon[M]^2$.
    Store the values of $\epsilon[M]^2$ in a vector |err_ldct|.
    Compares the Fourier, Wavelets, DCT and local-DCT approximations.
    """
    cR = sort(abs(fL(: )))
    if cR(n^2) >cR(1)
        cR = reverse(cR); % be sure it is in reverse order
    err_ldct = max(norm(f(: ))^2 - cumsum(cR.^2), 1e-10)
    clf; hold('on')
    h = plot(log10(err_fft / norm(f(: ))^2), 'r')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    h = plot(log10(err_dct / norm(f(: ))^2), 'g')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    h = plot(log10(err_ldct / norm(f(: ))^2), 'b')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    h = plot(log10(err_wav / norm(f(: ))^2), 'k')
    if using_matlab()  
        set(h, 'LineWidth', lw)
    axis(ax)
    title('log_{10}(\epsilon^2[M]/ ||f||^2)')
    legend('Fourier', 'DCT', 'local-DCT', 'Wavelets')
    if using_matlab()  
        box('on')


def exo12():
    """
    Compare the approximation error decay for those images.
    Display $ \log_{10}(\norm{f-f_M}) $ as a function of $\log_{10}(M)$.
    """
    err = []
    for i in 1: size(fList, 3):
        % 
        fW = perform_wavortho_transf(fList(: , : , i), Jmin, + 1, options)
        e = compute_l2err(fW, 1)
        err = [err e]
    sel = 10: round(n*n/ 10)
    Err = err(sel, : )
    Err = Err ./ repmat(Err(1, : ), [size(Err, 1), 1])
    h = plot(log10(sel), log10(Err))
    if using_matlab()  
        set(h, 'LineWidth', lw)
    title('log_{10}(\epsilon^2[M]/ ||f||^2)')
    legend('regular', 'phantom', 'lena', 'mandrill')
    axis([log10(sel(1)) log10(sel(length(sel))) -7 0])
    if using_matlab()  
        box('on')


