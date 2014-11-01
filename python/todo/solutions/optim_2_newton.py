def exo1():
    """
    Implement the Newton algorithm.
    Display the evolution of $f(x^{(\ell)})$ and $\norm{x^{(\ell)}-x^{(+\infty)}}$
    during the iterations.
    isplay
    """
    niter = 10
    x = [-1.5; 2.5]
    x = [1.7; 2.7]
    E = []; X = []
    for i in 1: niter:
        X(: , end + 1) = x
        x = x - pinv(Hessf(x))*Gradf(x)
        E(i) = f(x(1), x(2))
    myplot = lambda E: plot(log10(E(E >eps)))
    Xs = [1; 1]
    e = sqrt(sum((X-repmat(Xs, [1 niter])).^2))
    subplot(2, 1, 1)
    myplot(E);  axis tight
    title('log_{10}|E(x_k)-E^*|')
    subplot(2, 1, 2)
    myplot(e);  axis tight
    title('log_{10}|x_k-x^*|')


def exo2():
    """
    Display the evolution of $x^{(\ell)}$, from several starting points.
    """
    cols = {'r' 'g' 'b' 'c' 'm' 'y'}
    Xinit = {[-1.5; 2.5] [1.7; 2.7] [-.3; .85]}
    niter = 10
    clf; hold on
    imagesc(x2, x1, perform_hist_eq(F, 'linear'))
    colormap jet(256)
    for k in 1: length(Xinit):
        x = Xinit{k}; % [4*rand-2; 3.5*rand-.5]
        X = []
    for i in 1: niter:
            X(: , end + 1) = x
            x = x - pinv(Hessf(x))*Gradf(x)
        plot(X(2, : ), X(1, : ), [cols{k} '.-'], 'MarkerSize', 20, 'LineWidth', 2)
    axis([-.5 3 -2 2])


def exo3():
    """
    Implement the Newton descent algorithm.
     d = Hinv(Gradf(x), A(x), d);
    """
    niter = 12; % 300
    x = y
    E = []
    d = zeros(n)
    for i in 1: niter:
        E(end + 1) = f(x)
        [d, ~] = Hinv(Gradf(x), A(x), d)
        d = flatI(d)
        x = x - d
    h = plot(E)
    title('f(x^{l})')
    set(h, 'LineWidth', 2)
    axis tight


def exo4():
    """
    Compare the Newton descent with the gradient descent with a fixed step
    size, in term of decay of the energy.
    
    """
    tau = 1.8/ (1 + lambda*8/ epsilon)
    tau = tau*4
    x = y
    E1 = []
    for i in 1: niter:
        E1(end + 1) = f(x)
        x = x - tau*Gradf(x)
    h = plot([E; E1]')
    title('f(x^{l})')
    legend('Newton', 'Gradient')
    set(h, 'LineWidth', 2)
    axis tight


def exo5():
    """
    The direct comparison between gradient method and Newton is not fair in
    term of iteration count. Indeed, and iteration of Newton requires
    several steps of conjugate gradient, which takes some time.
    Try to set-up a fair comparison benchmark that takes into account the
    runing time of the methods. Pay a particular attention to the number of
    steps (or the tolerance criterion) that parameterize |cgs|.
    """


