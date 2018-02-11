def exo1():
    """
    Implement a full wavelet transform that extract iteratively wavelet
    coefficients, by repeating these steps. Take care of choosing the
    correct number of steps.
    """
    Jmin = 0
    fw = f
    for j in J: -1: Jmin:
        fw(1: 2^(j + 1), 1: 2^(j + 1)) = haar(fw(1: 2^(j + 1), 1: 2^(j + 1)))
        %
        j1 = J-j
        if j1 <4
            A = fw(1: 2^(j + 1), 1: 2^(j + 1))
            imageplot(A(1: 2^j, 2^j + 1: 2^(j + 1)), ['Horizontal, j = ' num2str(j)], 3, 4, j1 + 1)
            imageplot(A(2^j + 1: 2^(j + 1), 1: 2^j), ['Vertical, j = ' num2str(j)], 3, 4, j1 + 5)
            imageplot(A(2^j + 1: 2^(j + 1), 2^j + 1: 2^(j + 1)), ['Diagonal, j = ' num2str(j)], 3, 4, j1 + 9)


def exo2():
    """
    Write the inverse wavelet transform that computes $f_1$ from
    coefficients |fW|.
    """
    f1 = fw
    for j in Jmin: J:
        s = 1: 2^j; t = 2^j + 1: 2^(j + 1); u = 1: 2^(j + 1)
        f1(u, u) = ihaar(f1(s, s), f1(s, t), f1(t, s), f1(t, t))
        %
        j1 = J-j
        if j1 >0 & j1 <5
            A = f1(1: 2^(j + 1), 1: 2^(j + 1))
            subplot(2, 2, j1)
            imageplot(A, ['Partial reconstruction, j = ' num2str(j)])


def exo3():
    """
    Display the reconstructed signal obtained from |fw1|, for a decreasing cut-off scale $j$.
    """
    jlist = J-(1: 4)
    fw = perform_haar_transf(f, 1, + 1)
    for i in 1: length(jlist):
        j = jlist(i)
        fw1 = zeros(n); fw1(1: 2^j, 1: 2^j) = fw(1: 2^j, 1: 2^j)
        f1 = perform_haar_transf(fw1, 1, -1)
        % display
        subplot(2, 2, i)
        imageplot(f1)
        title(strcat(['j = ' num2str(j) ', SNR = ' num2str(snr(f, f1), 3) 'dB']))


def exo4():
    """
    Find the threshold $T$ so that the number of remaining coefficients in
    |fwT| is a fixed number $m$. Use this threshold to compute |fwT| and then display
    the corresponding approximation $f_1$ of $f$. Try for an increasing number $m$ of coeffiients.
    """
    m_list = round([.005 .01 .05 .1]*N); % number of kept coefficients
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
        subplot(2, 2, i)
        imageplot(f1)
        title(strcat(['m = ' num2str(m) ', SNR = ' num2str(snr(f, f1), 3) 'dB']))


