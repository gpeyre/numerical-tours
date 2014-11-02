def exo1():
    """
    Implement a full wavelet Haar transform that extract iteratively wavelet
    coefficients, by repeating these steps. Take care of choosing the
    correct number of steps.
    """
    J = log2(N)-1
    Jmin = 0
    fw = f
    subplot(4, 1, 1)
    plot(f); axis('tight'); title('Signal')
    for j in J: -1: Jmin:
        fw(1: 2^(j + 1)) = haar(fw(1: 2^(j + 1)))
        %    
        j1 = J-j
        if j1 <3
            d = fw(2^j + 1: 2^(j + 1))
            subplot(4, 1, j1 + 2)
            plot(1: 2^(j1 + 1): N, d);  axis('tight')
            title(strcat(['Details, j = ' num2str(j)]))


def exo2():
    """
    Write the inverse wavelet transform that computes |f1| from the
    coefficients |fw|.
    """
    f1 = fw
    for j in Jmin: J:
        f1(1: 2^(j + 1)) = ihaar(f1(1: 2^j), f1(2^j + 1: 2^(j + 1)))
        j1 = J-j
        if j1 <4
            subplot(4, 1, j1 + 1)
            plot(1: 2^j1: N, f1(1: 2^(j + 1)), '.-'); axis('tight')
            title(strcat(['Partial reconstruction, j = ' num2str(j)]))


def exo3():
    """
    Display the reconstructed signal obtained from |fw1|, for a decreasing cut-off scale $j$.
    """
    jlist = J-(1: 3)
    fw = perform_haar_transf(f, 1, + 1)
    for i in 1: length(jlist):
        j = jlist(i)
        fw1 = fw; fw1(2^j + 1: end) = 0
        f1 = perform_haar_transf(fw1, 1, -1)
        % display
        subplot(length(jlist), 1, i)
        hh = plot(f1); axis('tight')
        if using_matlab()
            set_linewidth(hh, 2)
        title(strcat(['j = ' num2str(j) ', SNR = ' num2str(snr(f, f1), 3) 'dB']))


def exo4():
    """
    Find the threshold $T$ so that the number of remaining coefficients in
    |fwT| is a fixed number $m$. Use this threshold to compute |fwT| and then display
    the corresponding approximation $f_1$ of $f$. Try for an increasing number $m$ of coeffiients.
    ompute the threshold T
    """
    m_list = round([.05 .1 .2]*N); % number of kept coefficients
    fw = perform_haar_transf(f, 1, + 1)
    for i in 1: length(m_list):
        m = m_list(i)
        % select threshold
        v = sort(abs(fw(: )))
        if v(1) <v(N)
            v = reverse(v)
        T = v(m)
        fwT = fw .* (abs(fw) >= T)
        % inverse
        f1 = perform_haar_transf(fwT, 1, -1)
        % display
        subplot(length(m_list), 1, i)
        hh = plot(f1); axis('tight')
        if using_matlab()
            set_linewidth(hh, 2)
        title(strcat(['m = ' num2str(m) ', SNR = ' num2str(snr(f, f1), 3) 'dB']))


def exo5():
    """
    Compute wavelets at several positions and scales.
    """
    J = log2(N)-1
    selj = (J-2: J)-3
    pos = [0 .5]
    f = []
    k = 0
    for j in selj:
        k = k + 1
    for q in 1: length(pos):
            fw = zeros(N, 1)
            p = 1 + (1 + pos(q))*2^j
            fw(p) = 1
            f(: , q) = perform_haar_transf(fw, 1, -1)
            f(: , q) = circshift(f(: , q), N/ 4)
        f(1: N/ 2-1, 2) = nan(); f(N/ 2 + 1: N, 1) = nan()
        subplot(3, 1, k)
        hh = plot(f); axis('tight')
        axis([1 N min(f(: ))*1.05 max(f(: ))*1.05])
        if using_matlab()
            set_linewidth(hh, 2)


