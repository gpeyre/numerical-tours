def exo1():
    """
    Perform the gradient descent using a fixed step size $\tau_k=\tau$.
    Display the decay of the energy $f(x^{(k)})$ through the iteration.
    Save the iterates so that |X(:,k)| corresponds to $x^{(k)}$.
    """
    x = [.5; .5]
    niter = 20
    E = []
    D = []
    X = []
    for i in 1: niter:
        X(: , i) = x
        E(end + 1) = f(x)
        D(end + 1) = norm(x)
        x = x - tau*Gradf(x)
    h = plot(log10(E))
    set(h, 'LineWidth', 2)
    axis tight
    title('log_{10}(x^{(k)})')


def exo2():
    """
    Display the iteration for several different step sizes.
    
    """
    taulist = [.3 1 1.7]/ eta
    xinit = [[.7; .7], [-.7; .5], [-.7; -.6]]
    collist = {'k' 'g' 'r'}
    clf; hold on
    imagesc(t, t, F); colormap jet(256)
    contour(t, t, F, 20, 'k')
    for k in 1: length(taulist):
        tau = taulist(k)
        %
        x = xinit(: , k)
        niter = 100
        X = []
    for i in 1: niter:
            X(: , i) = x
            x = x - tau*Gradf(x)
        %
        h = plot(X(1, : ), X(2, : ), [collist{k} '.-'])
        set(h, 'LineWidth', 2)
        set(h, 'MarkerSize', 15)
        axis off; axis equal


def exo3():
    """
    Implement the gradient descent. Monitor the decay of $f$ through the
    iterations.
    """
    niter = 300
    x = y
    E = []
    for i in 1: niter:
        E(end + 1) = f(y, x, epsilon)
        x = x - tau*Gradf(y, x, epsilon)
    h = plot(E)
    set(h, 'LineWidth', 2)
    axis tight


def exo4():
    """
    Display the evolution of the inpainting process.
    """
    tau = 1.8/ (8/ epsilon)
    tau = tau * 100
    E = []
    niter = 20000
    ndisp = round(linspace(1, niter, 5)); ndisp(1) = []
    x = y
    clf; q = 1
    for i in 1: niter:
        E(i) = J(x, epsilon)
        if i >1 && E(i) >E(i-1)
           tau = tau*.8
        x = x - tau * GradJ(x, epsilon)
        x = ProjH(x, y)
        if i = =ndisp(q)
            subplot(2, 2, q)
            imageplot(x)
            q = q + 1


def exo5():
    """
    Try with several values of $\epsilon$.
    au = tau * 100;
    """
    epsilon = 200
    tau = 1.8/ (8/ epsilon)
    E = []
    niter = 150
    ndisp = round(linspace(1, niter, 5)); ndisp(1) = []
    x = y
    clf; q = 1
    for i in 1: niter:
        E(i) = J(x, epsilon)
        if i >1 && E(i) >E(i-1)
           tau = tau*.8
        x = x - tau * GradJ(x, epsilon)
        x = ProjH(x, y)
        if i = =ndisp(q)
            subplot(2, 2, q)
            imageplot(x)
            q = q + 1


