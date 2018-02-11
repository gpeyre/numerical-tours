def exo1():
    """
    Compute an approximation |fLow| using the $ P=2^{2J}=(n/k_0)^2 $ low pass
    coefficients.
    """
    fwLow = zeros(n)
    fwLow(1: 2^J, 1: 2^J) = fw(1: 2^J, 1: 2^J)
    fLow = WavI(fwLow)
    myplot = lambda f1: imageplot(clamp(f1), ['PSNR = ' num2str(psnr(f, f1), 3) 'dB'])
    myplot(fLow)


def exo2():
    """
    Reconstruct an image using the pseudo inverse coefficients $\Phi^+ y =
    \Phi^* y$.
    """
    f1w = fw
    f1w(I0) = PhiS(y)
    fL2 = WavI(f1w)
    clf; myplot(fL2)


def exo3():
    """
    Implement the proximal and reversed-proximal mappings of $F$ (the orthogonal projector on
    $\Cc$ and $G$ (soft thresholding). In Matlab, use inline function with the |@|
    operator.
    """
    ProxF = lambda x, gamma: x + PhiS(y-Phi(x))
    ProxG = lambda x, gamma: max(0, 1-gamma./ max(1e-15, abs(x))).*x
    rProxF = lambda x, gamma: 2*ProxF(x, gamma)-x
    rProxG = lambda x, gamma: 2*ProxG(x, gamma)-x


def exo4():
    """
    Implement the DR iterative algorithm.
    Keep track of the evolution of the $\ell^1$ norm $G(x_k)$.
    """
    G = []
    F = []
    tx = zeros(N, 1)
    niter = 300
    for i in 1: niter:
        tx = (1-mu/ 2)*tx + mu/ 2*rProxG(rProxF(tx, gamma), gamma)
        x = ProxF(tx, gamma)
        G(i) = norm(x, 1)
        F(i) = norm(y-Phi(x))
    h = plot(G)
    set(h, 'LineWidth', 2)
    axis tight


def exo5():
    """
    Display the image reconstructed using the $P_0$ linear and $P$ CS
    measurements. The total number of used measurements is thus $P+P_0$.
    """
    fCSw = fw
    fCSw(I0) = x
    fCS = WavI(fCSw)
    clf; myplot(fCS)


def exo6():
    """
    define the proximal operator $ \text{Prox}_{\ga G} $ of $G$,
    and its reversed proximal mapping.
    """
    Energy = lambda x: sqrt(sum(x(I).^2))
    SoftAtten = lambda x, gamma: max(0, 1-gamma./ abs(x))
    EnergyAtten = lambda x, gamma: repmat(SoftAtten(Energy(x), gamma), [w*w 1])
    Flatten = lambda x: x(: )
    ProxG = lambda x, gamma: accumarray(I(: ), Flatten(EnergyAtten(x, gamma)), [N 1], @prod) .* x
    rProxG = lambda x, gamma: 2*ProxG(x, gamma)-x


def exo7():
    """
    Implement the DR iterative algorithm.
    Keep track of the evolution of $G(x_k)$.
    """
    g = []
    tx = zeros(N, 1)
    niter = 100
    for i in 1: niter:
        tx = (1-mu/ 2)*tx + mu/ 2*rProxG(rProxF(tx, gamma), gamma)
        x = ProxF(tx, gamma)
        g(i) = G(x)
    h = plot(g)
    set(h, 'LineWidth', 2)
    axis tight


def exo8():
    """
    Display the image reconstructed using the $P_0$ linear and $P$ CS
    measurements.
    """
    fCSBlockw = fw
    fCSBlockw(I0) = x
    fCSBlock = WavI(fCSBlockw)
    clf; myplot(fCSBlock)


