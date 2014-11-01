def exo1():
    """
    Perform the curve evolution.
    You need to resample it a few times.
    """
    gamma = gamma1
    displist = round(linspace(1, niter, 10))
    k = 1
    clf; hold on
    for i in 1: niter:
        gamma = resample(gamma + dt * normalC(gamma))
        if i = =displist(k)
            % display
            h = plot(gamma([1: end 1]), 'r')
            if i = =1 || i = =niter
                set(h, 'LineWidth', 2)
            k = k + 1
            drawnow
            axis('tight');  axis('off')


def exo2():
    """
    Perform the curve evolution.
    """
    gamma = gamma0
    displist = round(linspace(1, niter, 10))
    k = 1
    clf; hold on
    imageplot(W)
    for i in 1: niter:
        N = normal(gamma)
        g = EvalW(gamma).*normalC(gamma) - dotp(EvalG(gamma), N) .* N
        gamma = resample(gamma + dt*g)
        if i = =displist(k)       
            h = plot(imag(gamma([1: end 1])), real(gamma([1: end 1])), 'r')
            if i = =1 || i = =niter
                set(h, 'LineWidth', 2)
            k = k + 1
            drawnow
            axis('ij'); axis('off')


def exo3():
    """
    Create an initial circle $\gamma_0$ of $p$ points.
    """
    r = .95*n/ 2
    p = 128; % number of points on the curve
    theta = linspace(0, 2*pi, p + 1)'; theta(end) = []
    gamma0 = n/ 2*(1 + 1i) +  r*(cos(theta) + 1i*sin(theta))
    gamma = gamma0
    clf; hold on
    imageplot(f)
    h = plot(imag(gamma([1: end 1])), real(gamma([1: end 1])), 'r')
    set(h, 'LineWidth', 2)


def exo4():
    """
    Perform the curve evolution.
    
    """
    options.order = 2
    G = grad(W, options)
    G = G(: , : , 1) + 1i*G(: , : , 2)
    EvalG = lambda gamma: interp2(1: n, 1: n, G, imag(gamma), real(gamma))
    EvalW = lambda gamma: interp2(1: n, 1: n, W, imag(gamma), real(gamma))
    gamma = gamma0
    displist = round(linspace(1, niter, 10))
    k = 1
    clf; hold on
    imageplot(f)
    for i in 1: niter:
        n = normal(gamma)
        g = EvalW(gamma).*normalC(gamma) - dotp(EvalG(gamma), n) .* n
        gamma = resample(gamma + dt*g)
        if i = =displist(k)       
            h = plot(imag(gamma([1: end 1])), real(gamma([1: end 1])), 'r')
            if i = =1 || i = =niter
                set(h, 'LineWidth', 2)
            k = k + 1
            drawnow
            axis('ij'); axis('off')


def exo5():
    """
    Compute an edge attracting criterion $W(x)>0$, that is small in area of strong
    gradient.
    """
    options.order = 2
    G = grad(f, options)
    G = sqrt(sum(G.^2, 3))
    G = perform_blurring(G, 3)
    G = min(G, .4)
    W = rescale(-G, .4, 1)
    imageplot(W)


def exo6():
    """
    Perform the curve evolution.
    Be careful to impose the boundary conditions at each step.
    
    """
    options.order = 2
    G = grad(W, options)
    G = G(: , : , 1) + 1i*G(: , : , 2)
    EvalG = lambda gamma: interp2(1: n, 1: n, G, imag(gamma), real(gamma))
    EvalW = lambda gamma: interp2(1: n, 1: n, W, imag(gamma), real(gamma))
    gamma = gamma0
    displist = round(linspace(1, niter, 10))
    k = 1
    clf; hold on
    imageplot(f)
    for i in 1: niter:
        N = normal(gamma)
        g = EvalW(gamma).*normalC(gamma) - dotp(EvalG(gamma), N) .* N
        gamma = resample(gamma + dt*g)
        % impose start/ end point
        gamma(1) = x0; gamma(end) = x1
        if i = =displist(k)   
            h = plot(imag(gamma([1: end])), real(gamma([1: end])), 'r')
            if i = =1 || i = =niter
                set(h, 'LineWidth', 2)
            h = plot(imag(gamma([1 end])), real(gamma([1 end])), 'b.'); set(h, 'MarkerSize', 30)
            axis('ij'); axis('off')
            k = k + 1
            drawnow


