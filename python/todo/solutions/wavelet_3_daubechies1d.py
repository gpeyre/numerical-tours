def exo1():
    """
    Implement a full wavelet transform that extract iteratively wavelet
    coefficients, by repeating these steps. Take care of choosing the
    correct number of steps.
    """
    Jmax = log2(N)-1; Jmin = 0
    fw = f
    subplot(4, 1, 1)
    plot(f); axis('tight'); title('Signal')
    for j in Jmax: -1: Jmin:
        a1 = fw(1: 2^(j + 1))
        a = subsampling(cconvol(a1, h))
        d = subsampling(cconvol(a1, g))
        fw(1: 2^(j + 1)) = cat(1, a, d)
        j1 = Jmax-j
        if j1 <3
            subplot(4, 1, j1 + 2)
            plot(1: 2^(j1 + 1): N, d);  axis('tight')
            title(['Details, j = ' num2str(j)])


def exo2():
    """
    Write the inverse wavelet transform that computes |f1| from the
    coefficients |fw|.
    """
    f1 = fw
    for j in Jmin: Jmax:
        a = f1(1: 2^j)
        d = f1(2^j + 1: 2^(j + 1))
        a = cconvol(upsampling(a, 1), reverse(h), 1)
        d = cconvol(upsampling(d, 1), reverse(g), 1)
        f1(1: 2^(j + 1)) = a + d
        j1 = Jmax-j
        if j1 <4
            subplot(4, 1, j1 + 1)
            plot(1: 2^j1: N, f1(1: 2^(j + 1)), '.-'); axis('tight')
            title(['Partial reconstruction, j = ' num2str(j)])


def exo3():
    """
    Find the threshold $T$ to obtained a given number $M$ of
    non thresholded coefficients.
    Try for an increasing number $M$ of coeffiients.
    ompute the threshold T
    """
    f = load_signal(name, N)
    m_list = [25 50 100]; % number of kept coefficients
    options.h = compute_wavelet_filter('Daubechies', 4)
    for i in 1: length(m_list):
        m = m_list(i)
        fw = perform_wavortho_transf(f, 1, + 1, options)
        % select threshold
        v = sort(abs(fw(: )))
        if v(1) <v(N)
            v = reverse(v)
        T = v(m)
        fwT = fw .* (abs(fw) >= T)
        % inverse
        f1 = perform_wavortho_transf(fwT, 1, -1, options)
        % display
        subplot(length(m_list), 1, i)
        hh = plot(f1); axis('tight')
        if using_matlab()
            set_linewidth(hh, 2)
        title(['M = ' num2str(m) ', SNR = ' num2str(snr(f, f1), 3) 'dB'])


def exo4():
    """
    Try with
    Different kind of wavelets, with an increasing number of vanishing
    moments.
    ompute the threshold T
    """
    m = 100; % number of kept coefficients
    vmlist = 1: 3
    for i in 1: length(vmlist):
        vm = vmlist(i)
        options.h = compute_wavelet_filter('Daubechies', vm*2)
        fw = perform_wavortho_transf(f, 1, + 1, options)
        % select threshold
        v = sort(abs(fw(: )))
        if v(1) <v(N)
            v = reverse(v)
        T = v(m)
        fwT = fw .* (abs(fw) >= T)
        % inverse
        f1 = perform_wavortho_transf(fwT, 1, -1, options)
        % display
        subplot(length(vmlist), 1, i)
        hh = plot(f1); axis('tight')
        if using_matlab()
            set_linewidth(hh, 2)
        title([num2str(vm) ' VM, SNR = ' num2str(snr(f, f1), 3) 'dB'])


def exo5():
    """
    Compute wavelets at several positions and scales.
    """
    Jmax = log2(N)-1
    options.h = compute_wavelet_filter('Daubechies', 2*2)
    selj = (Jmax-2: Jmax)-3
    pos = [0 .5]
    f = []
    k = 0
    for j in selj:
        k = k + 1
    for q in 1: length(pos):
            fw = zeros(N, 1)
            p = 1 + (1 + pos(q))*2^j
            fw(p) = 1
            f(: , q) = perform_wavortho_transf(fw, 1, -1, options)
            f(: , q) = circshift(f(: , q), N/ 4)
        f(1: N/ 2-1, 2) = nan(); f(N/ 2 + 1: N, 1) = nan()
        subplot(3, 1, k)
        hh = plot(f); axis('tight')
        axis([1 N min(f(: ))*1.05 max(f(: ))*1.05])
        if using_matlab()
            set_linewidth(hh, 2)


def exo6():
    """
    Display Daubechies wavelets with an increasing number of vanishing
    moments.
    """
    vmlist = 2: 5
    j = 3
    for i in 1: length(vmlist):
        vm = vmlist(i)
        options.h = compute_wavelet_filter('Daubechies', vm*2)
        fw = zeros(N, 1)
        fw(1 + 2^j) = 1
        f = perform_wavortho_transf(fw, 1, -1, options)
        f = circshift(f, N/ 2)
        subplot(length(vmlist), 1, i)
        hh = plot(f); axis('tight')
        if using_matlab()
            set_linewidth(hh, 2)
        title([num2str(vm) ' VM'])


