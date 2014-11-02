def exo1():
    """
    Compute and display $d_y$
    VD
    _y
    isp
    """
    [U, S, V] = svd(MusicHankel(y), 0); S = diag(S)
    Ubot = U(: , N + 1: end)
    d = Ubot'*exp(-2i*pi*(0: L-1)'*z(: )')
    d = sum(abs(d).^2) / L
    clf; hold on
    plot(z, d, 'b')
    stem(x0, 1 + x0*0, 'r.')
    axis([0 1 0 1]); box on


def exo2():
    """
    Display the roots of $P_y$ that are inside the unit disk
    oefficients
    oots
    eep those inside
    isplay
    """
    B = []
    for j in 1: size(Ubot, 2):
        u = Ubot(: , j)
        v = flipud(conj(u))
        B(: , j) = conv(u, v)
    C = sum(B, 2)
    R = roots(C(end: -1: 1))
    R = R(abs(R) <= 1)
    clf; hold on
    plot(exp(2i*pi*z), 'k')
    plot(R, 'b.', 'MarkerSize', ms)
    axis equal; box on
    axis([-1 1 -1 1]*1.1)
    axis off


def exo3():
    """
    Keep only the best $N$ ones, i.e. those that are the closest from the
    unit circle.
    We denote those as $\tilde x \in \mathbb{T}^N$.
    Recover an approximation $\tilde a$ of the amplitudes $a_0$
    By solving in least squares (using backslash operator) the system
    $ \Phi_{\tilde x} \tilde a = y. $
    Display the recovered measure $ m_{\tilde a, \tilde x} $.
    osition
    eep only the best N ones.
    ompute amplitude by solving a least square.
    isplay the recovered measure.
    """
    [x1, I] = sort(mod(angle(R), 2*pi)/ (2*pi))
    R = R(I)
    [~, I] = sort(abs(abs(R)-1))
    x1 = x1(I(1: N))
    R = R(I(1: N))
    a1 = real(Phi(x1)\y)
    clf; hold on
    plot(z, f)
    stem(x0, a0, 'r.')
    stem(x1, a1, 'k--')


def exo4():
    """
    Display the evolution of roots as the noise level $\sigma$ increases.
    
    
    
    """
    slist = linspace(1e-9, 1, 1000)
    RL = []
    for is in 1: length(slist):
        % observations
        sigma = slist(is)
        y = y0 + sigma*norm(y0)*w
        % SVD
        [U, S, V] = svd(MusicHankel(y), 0); S = diag(S)
        Ubot = U(: , N + 1: end)
        % Coefficients
        B = []
    for j in 1: size(Ubot, 2):
            u = Ubot(: , j)
            v = flipud(conj(u))
            B(: , j) = conv(u, v)
        C = sum(B, 2)
        % Roots
        R = roots(C(end: -1: 1))
        % keep those inside
        R = R(abs(R) <= 1)
        % position
        [x1, I] = sort(mod(angle(R), 2*pi)/ (2*pi))
        R = R(I)
        % Keep only the best N ones.
        if 0
            [~, I] = sort(abs(abs(R)-1))
            x1 = x1(I(1: N))
        RL(: , end + 1) = R
    clf; hold on
    plot(exp(2i*pi*z), 'k')
    cm = jet(length(slist))
    for is in 1: length(slist):
        plot(RL(: , is), '.',  'MarkerSize', ms, 'color', cm(is, : ))
    axis equal; box on
    axis([-1 1 -1 1]*1.1)
    axis off


