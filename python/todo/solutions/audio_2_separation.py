def exo1():
    """
    Compute the STFT of the micros, and store them into a matrix |Y|.
    """
    for i in 1: p:
        Y(: , : , i) = perform_stft(y(: , i), w, q, options)
        subplot(p, 1, i)
        plot_spectrogram(Y(: , : , i))
        set_graphic_sizes([], 20)
        title(strcat('Micro #', num2str(i)))


def exo2():
    """
    Display some points of |P| in the transformed (time/frequency) domain.
    """
    sel = randperm(size(P, 1)); sel = sel(1: npts)
    plot(P(sel, 1), P(sel, 2), '.')
    axis([-1 1 -1 1]*5)
    set_graphic_sizes([], 20)
    title('Transformed domain')


def exo3():
    """
    The histogram computed from the whole set of points are not peacked
    enough. To stabilize the detection of mixing direction, compute an
    histogram from a reduced set of point that have the largest amplitude.
    ompute the energy of each point
    xtract only a small sub-set
    ompute histogram
    isplay histograms
    """
    d = sum(P.^2, 2)
    rho = .1
    [v, I] = sort(d)
    if v(1) <v(length(v))
        I = reverse(I)
    I = I(1: round(rho*length(I)))
    P1 = P(I, : );  % transformed points
    Theta = mod(atan2(P1(: , 2), P1(: , 1)), pi())
    nbins = 200
    [h, t] = hist(Theta, nbins)
    h = h/ sum(h)
    bar(t, h); axis('tight')


def exo4():
    """
    Detect the direction |M1| approximating the true direction |M| by
    looking at the local maxima of the histogram. First detect the set of
    local maxima, and then keep only the three largest.
    ort in descending order
    """
    s1 = [2: nbins nbins-1]
    s2 = [2 1: nbins-1]
    I = find((h(s1) <h) & (h(s2) <h))
    [v, u] = sort(h(I))
    if v(1) <v(length(v))
        u = reverse(u)
    theta1 = t(I(u(1: 3))); theta1 = theta1(: )'
    M1 = [cos(theta1); sin(theta1)]
    M, M1


