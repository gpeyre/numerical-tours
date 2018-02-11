def exo1():
    """
    Implement a full wavelet transform that extract iteratively wavelet
    coefficients, by repeating these steps. Take care of choosing the
    correct number of steps.
    """
    Jmax = log2(n)-1
    Jmin = 0
    fW = f
    for j in Jmax: -1: Jmin:
        A = fW(1: 2^(j + 1), 1: 2^(j + 1))
    for d in 1: 2:
            Coarse = subsampling(cconvol(A, h, d), d)
            Detail = subsampling(cconvol(A, g, d), d)
            A = cat3(d, Coarse, Detail)
        fW(1: 2^(j + 1), 1: 2^(j + 1)) = A
        j1 = Jmax-j
        if j1 <4
            imageplot(A(1: 2^j, 2^j + 1: 2^(j + 1)), strcat(['Horizontal, j = ' num2str(j)]), 3, 4, j1 + 1)
            imageplot(A(2^j + 1: 2^(j + 1), 1: 2^j), strcat(['Vertical, j = ' num2str(j)]), 3, 4, j1 + 5)
            imageplot(A(2^j + 1: 2^(j + 1), 2^j + 1: 2^(j + 1)), strcat(['Diagonal, j = ' num2str(j)]), 3, 4, j1 + 9)


def exo2():
    """
    Write the inverse wavelet transform that computes |f1| from the
    coefficients |fW|. Compare |f1| with |f|.
    """
    f1 = fW
    for j in Jmin: Jmax:
        A = f1(1: 2^(j + 1), 1: 2^(j + 1))
    for d in 1: 2:
            if d = =1
                Coarse = A(1: 2^j, : )
                Detail = A(2^j + 1: 2^(j + 1), : )
            else
                Coarse = A(: , 1: 2^j)
                Detail = A(: , 2^j + 1: 2^(j + 1))
            Coarse = cconvol(upsampling(Coarse, d), reverse(h), d)
            Detail = cconvol(upsampling(Detail, d), reverse(g), d)
            A = Coarse + Detail
            j1 = Jmax-j
            if j1 >0 & j1 <5
                subplot(2, 2, j1)
                imageplot(A, strcat(['Partial reconstruction, j = ' num2str(j)]))
        f1(1: 2^(j + 1), 1: 2^(j + 1)) = A


def exo3():
    """
    Compute and display the linear approximation |fLin| obtained from the
    coefficients |fWLin| by inverse wavelet transform.
    """
    fLin = perform_wavortho_transf(fWLin, 0, -1)
    elin = snr(f, fLin)
    imageplot(f, 'Original', 1, 2, 1)
    imageplot(clamp(fLin), strcat(['Linear, SNR = ' num2str(elin)]), 1, 2, 2)


def exo4():
    """
    Find the thresholds |T| so that the number |m| of remaining coefficients in
    |fWT| are |m=n^2/16|. Use this threshold to compute |fWT| and then display
    the corresponding approximation |Mnlin| of |f|. Compare this result with
    the linear approximation.
    umber of kept coefficients
    ompute the threshold T
    elect threshold
    nverse transform
    nverse
    isplay
    """
    m = round(n^2/ 16)
    Jmin = 1
    fW = perform_wavortho_transf(f, Jmin, + 1)
    v = sort(abs(fW(: )))
    if v(1) <v(n^2)
        v = reverse(v)
    T = v(m)
    fWT = fW .* (abs(fW) >= T)
    Mnlin = perform_wavortho_transf(fWT, Jmin, -1)
    imageplot(clamp(fLin), strcat(['Linear, SNR = ' num2str(snr(f, fLin), 3) 'dB']), 1, 2, 1)
    imageplot(clamp(Mnlin), strcat(['Non-linear, SNR = ' num2str(snr(f, Mnlin), 3) 'dB']), 1, 2, 2)


def exo5():
    """
    Try with
    Different kind of wavelets, with an increasing number of vanishing
    moments.
    """
    vmlist = 1: 3
    for i in 1: length(vmlist):
        vm = vmlist(i)
        options.h = compute_wavelet_filter('Daubechies', vm*2)
        fW = perform_wavortho_transf(f, 1, + 1, options)
        % select threshold
        v = sort(abs(fW(: )))
        if v(1) <v(n)
            v = reverse(v)
        T = v(m)
        fWT = fW .* (abs(fW) >= T)
        % inverse
        f1 = perform_wavortho_transf(fWT, 1, -1, options)
        % display
        imageplot(clamp(f1), strcat([num2str(vm) ' VM, SNR = ' num2str(snr(f, f1), 3) 'dB']), 1, 3, i)


def exo6():
    """
    Implement the foward separable transform.
    Wavelet transformm in 1D each column |f(:,i)| to obtain coefficients |fWSep(:,i)|.
    Then re-transform each row |fWSep(i,:)'|, and store the result in |fW(i,:)'|.
    """
    fWSep = f
    Jmin = 1
    for i in 1: n:
        fWSep(: , i) = perform_wavortho_transf(fWSep(: , i), Jmin, + 1)
    for i in 1: n:
        fWSep(i, : ) = perform_wavortho_transf(fWSep(i, : )', Jmin, + 1)'


def exo7():
    """
    Implement the backward separable transform to recover an image |f1|
    from the coefficients |fWSep|, which backward transform each row and then each columns.
    """
    f1 = fWSep
    for i in 1: n:
        f1(: , i) = perform_wavortho_transf(f1(: , i), Jmin, -1)
    for i in 1: n:
        f1(i, : ) = perform_wavortho_transf(f1(i, : )', Jmin, -1)'


def exo8():
    """
    Perform |m=n^2/16|-terms best approximation with separable wavelet transform.
    Compare the result with isotropic wavelet approximation.
    isplay
    """
    f1 = perform_thresholding(fWSep, m, 'largest')
    for i in 1: n:
        f1(: , i) = perform_wavortho_transf(f1(: , i), Jmin, -1)
    for i in 1: n:
        f1(i, : ) = perform_wavortho_transf(f1(i, : )', Jmin, -1)'
    imageplot(clamp(Mnlin), strcat(['Isotropic, SNR = ' num2str(snr(f, Mnlin), 3) 'dB']), 1, 2, 1)
    imageplot(clamp(f1), strcat(['Separable, SNR = ' num2str(snr(f, f1), 3) 'dB']), 1, 2, 2)


def exo9():
    """
    Compute wavelets at several scales and orientation.
    Here we show only horizontal wavelets, in 2-D.
    """
    Jmax = log2(n)-1
    options.h = compute_wavelet_filter('Daubechies', 6)
    selj = (Jmax-2: Jmax)-3
    k = 0
    for j in selj:
    for s in 1: 3:
            k = k + 1
            fW = zeros(n, n)
            MWj = zeros(2^j, 2^j)
            MWj(1, 1) = 1
            if s = =1
                fW(1: 2^j, 2^j + 1: 2*2^j) = MWj
                str = 'Hor'
            elseif s = =2
                fW(2^j + 1: 2*2^j, 1: 2^j) = MWj
                str = 'Vert'
            elseif s = =3
                fW(2^j + 1: 2*2^j, 2^j + 1: 2*2^j) = MWj
                str = 'Diag'
            f1 = perform_wavortho_transf(fW, 1, -1, options)
            f1 = circshift(f1, [n n]/ 2)
            imageplot(f1, strcat(['j = ' num2str(j) ' ' str]), 3, 3, k)


def exo10():
    """
    Display Daubechies wavelets with different orientation, for different number of VM.
    """
    vmlist = [2 5]
    j = 3
    k = 0
    for i in 1: length(vmlist):
        vm = vmlist(i)
        options.h = compute_wavelet_filter('Daubechies', vm*2)
    for s in 1: 3:
            k = k + 1
            fW = zeros(n, n)
            if s = =1
                fW(1 + 2^j, 1) = 1
                str = 'Hor'
            elseif s = =2
                fW(1, 1 + 2^j) = 1
                str = 'Vert'
            elseif s = =3
                fW(1 + 2^j, 1 + 2^j) = 1
                str = 'Diag'
            f1 = perform_wavortho_transf(fW, 1, -1, options)
            [tmp, u] = max(f1(: )); [u, v] = ind2sub([n n], u)
            f1 = circshift(f1, [n/ 2-u n/ 2-v])
            subplot(2, 3, k)
            surf(f1)
            shading('interp'); view(3)
            axis('tight'); axis('off')
            lighting('phong')
            camlight
            colormap(jet(256))
            title(strcat([num2str(vm) ' VM, ' str]))


